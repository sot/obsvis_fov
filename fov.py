#!/usr/bin/env python

"""
========
fov
========

This code queries the sybase OCAT table (and as needed an ACA obspar table) to 
create obsvis-compatible field-of-view files for requested obsids, sequence
numbers, or proposals.
"""


import re
import os
import getpass
import urllib
import logging
from BeautifulSoup import BeautifulSoup

import asciitable
import numpy as np
import Ska.DBI
from Chandra.Time import DateTime


def default_offsets(det):
    """Get default offsets.

    :param det: one of HRC-I, HRC-S, ACIS-I, or ACIS-S
    :return: (y_offset, z_offset)
    :rtype: tuple
    """

    defaults = {'HRC-I': (0.0, 0.0),
                'HRC-S': (0.0, 0.0),
                'ACIS-I': (-0.2, -0.25),
                'ACIS-S': (0.15, -0.25)}
    return defaults[det]

last_aimpoint_shift = DateTime('2011:187')
default_lts = 'http://cxc.harvard.edu/target_lists/longsched.html'
logger = logging.getLogger()
console = logging.StreamHandler()
logger.setLevel(logging.INFO)
logger.addHandler(console)


def lts_table(lts_url=default_lts):
    """
    Grab the HTML long term schedule, skip everything except lines that
    look like obsid lines, and cast them into a table using asciitable.

    This does not get the lts week information available on that page from
    the paragraph order, but that can be retrieved directly from the
    axafocat.targets table

    :param lts_url: long term schedule url
    :returns: long term schedule targets in page order
    :rtype: numpy array

    """

    soup = BeautifulSoup(urllib.urlopen(lts_url).read())
    page = "".join(soup.findAll(text=True))
    pagelines = page.splitlines()
    obslines = []
    for line in pagelines:
        # the lines that begin with a sequence number and obsid and end with
        # proposal number
        if re.match('^\d{6}\s+(P1)?\s+\d{1,5}.*\d{8}\s*$', line):
            obslines.append(line)

    col_starts = (0, 7, 10, 16, 38, 45, 53, 61, 67,
                73, 80, 82, 84, 89, 110, 115,
                119, 122, 125, 128, 131, 134,
                138, 140, 165)
    col_ends = (6, 9, 14, 36, 43, 52, 60, 66, 72,
                79, 81, 83, 88, 107, 112, 117,
                120, 123, 126, 129, 132, 135,
                139, 162, 173)
    names = ('seqnbr', 'pool', 'obsid', 'target', 'ksec',
             'ra', 'dec', 'roll', 'pitch',
             'si', 'R', 'O', 'grating', 'observer', 'type', 'AO', 'OR_num',
             'TC', 'RC', 'PC', 'UC', 'Mlt',
             'CRem', 'approved', 'propid')

    lts = asciitable.read(obslines,
                          Reader=asciitable.FixedWidthNoHeader,
                          col_starts=col_starts,
                          col_ends=col_ends,
                          names=names
                      )
    return lts


def prelim_roll(obsid, lts=None):
    """Grab the roll from the long term schedule for the specified obsid."""
    if lts is None:
        lts = lts_table()
    match = lts['obsid'] == obsid
    n_match = len(np.flatnonzero(match))
    if n_match == 0:
        return None
    if n_match > 1:
        raise ValueError
    return lts[match]['roll'][0]


def archive_roll(obsid, dbh):
    """Grab the roll from the observations table for the specified obsid."""
    obs = dbh.fetchall("""select obsid, obi, ra_targ, dec_targ,
                          ra_pnt, dec_pnt, roll_pnt
                          from observations where obsid = %d""" % obsid)
    if not len(obs):
        return None
    if np.unique(obs['obi']) > 1:
        raise ValueError("Cannot determine roll for multi-obi obsid")
    if np.std(obs['roll_pnt'] > 2):
        raise ValueError("Observed rolls for obsid have std > 2.00; failing")
    # if there is more than one and they are close, just return the mean
    return np.mean(obs['roll_pnt'])


def get_obsids_seqnum(reqid, dbh):
    """Get the observations for a sequence from the target table."""
    targs = dbh.fetchall("select obsid from target where seq_nbr = '%d'"
                         % reqid)
    if len(targs):
        return targs['obsid']


def get_obsids_propnum(reqid, dbh):
    """Get the observations for a proposal from the ocat."""
    targs = dbh.fetchall("""select obsid from target t
join prop_info p on t.ocat_propid = p.ocat_propid
where prop_num = '%d'""" % reqid)
    if len(targs):
        return targs['obsid']


def get_fov(obsid, ocat_dbh, aca_dbh):
    """
    For a given obsid, retrieve the pieces necessary to build an obsvis FOV.

    :param obsid: requested obsid
    :param ocat_dbh: handle for axafocat database
    :param aca_dbh: handle for aca database
    :returns: field-of-view
    :rtype: dict
    """

    target = ocat_dbh.fetchone("select * from target where obsid = %d" % obsid)

    if not target:
        logger.warn("No entry for %d in target table" % obsid)
        return None

    roll = archive_roll(obsid, aca_dbh)
    if roll is None:
        lts = lts_table()
        roll = prelim_roll(obsid, lts)
    if roll is None:
        roll = 0

    # fov dictionary with mix of defaults and minimum info from target table
    fov = dict(fovid="obs%05d" % target['obsid'],
               groupid='',
               grouping='None',
               roll=roll,
               coordinates="%(ra)f, %(dec)f" % target,
               target=target['targname'],
               offsetsimz=0.0,
               offsetsimzmm=0.0,
               grating=target['grating'],
               offsety=target['y_det_offset'],
               offsetz=target['z_det_offset'],
               aca=0,
               instrument=target['instrument'])

    # many offsets seem to be NULL in the table, which is annoying.
    # this uses defaults for current display
    # as these have changed over time, this are not good for viewing
    # old observations
    if (fov['offsety'] is None) or (fov['offsetz'] is None):
        def_offsets = default_offsets(target['instrument'])
        if fov['offsety'] is None:
            fov['offsety'] = def_offsets[0]
        if fov['offsetz'] is None:
            fov['offsetz'] = def_offsets[1]

    sim = ocat_dbh.fetchone("select * from sim where obsid = %d" % obsid)
    if sim:
        arcmin_mm = 2.925
        # just set one or the other (simzmm or simz) or obsvis complains
        # fov['offsetsimzmm'] = sim['trans_offset']
        fov['offsetsimz'] = sim['trans_offset'] / arcmin_mm

    if (target['acisid']):
        fov['unselectedchips'] = 0
        fov['subarrays'] = None
        acis = ocat_dbh.fetchone("""select * from acisparam a, target t
                               where a.acisid=t.acisid
                               and a.acisid=%d""" % target['acisid'])
        i_ccd_map = {'ccdi%d_on' % i: 'I%d' % i for i in range(0, 4)}
        s_ccd_map = {'ccds%d_on' % i: 'S%d' % i for i in range(0, 6)}
        ccd_map = dict(i_ccd_map.items() + s_ccd_map.items())
        ccd_list = []
        for x in ccd_map:
            if acis[x] == 'Y' or re.match('O.+', acis[x]):
                ccd_list.append('%s' % ccd_map[x])
        fov['chips'] = ' '.join(ccd_list)
        if acis['subarray'] != 'NONE':
            fov['subarrays'] = \
"""Custom
      start row: %(subarray_start_row)s
      number of rows: %(subarray_row_count)s""" % acis
        fov['exposuremode'] = acis['exp_mode']
        fov['nodeboundaries'] = 0

    if (target['hrcid']):
        fov['blankingmode'] = 'None'
        fov['chipboundary'] = 1
        hrc = ocat_dbh.fetchone("""select * from hrcparam h, target t
                                 where h.hrcid=t.hrcid
                                 and h.hrcid=%d""" % target['hrcid'])
        if target['instrument'] == 'HRC-S':
            fov['unselectedchips'] = 0
            fov['timingmode'] = 0
            if hrc['timing_mode'] == 'Y':
                fov['timingmode'] = 1

    ot_date = None
    if target['lts_lt_plan'] is not None:
        ot = target['lts_lt_plan']
        ot_date = DateTime('%04d:%03d' % (ot.year, ot.dayofyear))

    if target['soe_st_sched_date'] is not None:
        ot = target['soe_st_sched_date']
        ot_date = DateTime('%04d:%03d' % (ot.year, ot.dayofyear))

    if ot_date and ot_date.secs < last_aimpoint_shift.secs:
        logger.warn("inaccurate fov: obsid %d predates last aimpoint shift."
                    % obsid)

    return fov


def fov_text(fov):
    """Given an fov dictionary return an array of lines of a fov file."""

    cols = ['groupid',
            'grouping',
            'coordinates',
            'target',
            'roll',
            'offsety',
            'offsetz',
            'offsetsimz',
            'offsetsimzmm',
            'grating',
            'aca',
            'instrument',
            'unselectedchips',
            'chips',
            'subarrays',
            'exposuremode',
            'nodeboundaries',
            'blankingmode',
            'chipboundary',
            'timingmode']

    fov_lines = []
    fov_lines.append("fovid: %(fovid)s" % fov)
    for col in cols:
        if col in fov:
            fov_lines.append("   %s: %s" % (col, fov[col]))
    return fov_lines


def main(reqids=[], ocat_user=None, password_file=None, outdir='.'):
    """
    For a list of ids, build field of views and write out to files.

    :param reqids: list of ids (obsids, sequence numbers, proposal ids)
    :param outdir: output directory name
    """
    # if both not defined, use defaults
    if ocat_user is None and password_file is None:
        username = getpass.getuser()
        dbpassfile = os.path.join(os.environ['HOME'], '.arc5gl_pwd')
    else:
        # if one or other defined, don't assume .arc5gl_pwd
        username = ocat_user
        if username is None:
            username = getpass.getuser()
        dbpassfile = password_file

    if dbpassfile and os.path.exists(dbpassfile):
        passwd = open(dbpassfile).read().strip()
    else:
        prompt = "No password-file found, supply password: "
        passwd = getpass.getpass(prompt=prompt)

    ocat_dbh = Ska.DBI.DBI(dbi='sybase', server='sqlsao',
                           user=username, passwd=passwd,
                           database='axafocat')

    aca_dbh = Ska.DBI.DBI(dbi='sybase', server='sybase',
                         user='aca_read', database='aca')

    obsids = []
    for reqid in reqids:
        len_reqid = len(str(reqid))
        if not ((len_reqid <= 5) or (len_reqid == 6) or (len_reqid == 8)):
            raise ValueError(
                "Requested id neither obsid, sequence number, nor proposal.")
        if len(str(reqid)) <= 5:
            logger.info("building fovs for obsid %d" % reqid)
            obsids.extend([reqid])
        if len(str(reqid)) == 6:
            logger.info("building fovs for obsids for seqnum %d" % reqid)
            sobs = get_obsids_seqnum(reqid, ocat_dbh)
            if sobs is not None:
                logger.info("\tobsids: %s" % ' '.join([str(x) for x in sobs]))
                obsids.extend(sobs)
            else:
                logger.info("\tNo obsids found for sequence")
        if len(str(reqid)) == 8:
            logger.info("building fovs for obsids for proposal %d" % reqid)
            sobs = get_obsids_propnum(reqid, ocat_dbh)
            if sobs is not None:
                logger.info("\tobsids: %s" % ' '.join([str(x) for x in sobs]))
                obsids.extend(sobs)
            else:
                logger.info("\tNo obsids found for proposal")

    for obsid in obsids:
        fov = get_fov(obsid, ocat_dbh, aca_dbh)

        if not fov:
            continue

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        outfile = os.path.join(outdir, 'obs%05d.fov' % obsid)
        fovfile = open(outfile, 'w')
        fov_lines = fov_text(fov)
        fovfile.writelines('\n'.join(fov_lines))
        fovfile.close()

        logger.info("obsid %d fov written to: %s" % (obsid, outfile))


def get_args():
    from argparse import ArgumentParser
    description = """
Obsvis FOV generator.  This tool generates obsvis-compatible field of view 
files.
It accepts as arguments a space-delimited list of obsids, sequence numbers,
or proposal ids (the argument type is determined by the number of digits).
For sequence numbers or proposal ids, more than one obsid may be retrieved,
and thus more than one field-of-view file will be created.

To create the field-of-view files, the targets table of the OCAT database is
queried.  Options are provided to assist with supplying a username/password
combination for access.  By default, the username of the user running the tool
(guessed from environment variables or uid) and the password in ~/.arc5gl_pwd
will be used.
"""

    parser = ArgumentParser(description=description)
    parser.add_argument('reqids',
                        type=int,
                        nargs='+',
                        help="id or ids to fetch")
    parser.add_argument("--ocat-user",
                        help="user for axafocat database (default=current user)")
    parser.add_argument("--password-file",
                        help="ocat user password file (default=~/.arc5gl_pwd)")
    parser.add_argument("--outdir",
                        default=".",
                        help="output directory for fov files")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    my_args = get_args()
    main(**my_args.__dict__)
