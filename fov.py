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
import logging
import Ska.Sun
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


def archive_roll(obsid, dbh):
    """Grab the roll from the obidet table for the specified obsid."""
    obs = dbh.fetchall("""select obsid, obi, roll
                          from axafapstat..obidet_0_5 where obsid = %d""" % obsid)
    if not len(obs):
        return None
    if np.std(obs['roll'] > 2):
        raise ValueError("Observed rolls for obsid have std > 2.00; failing")
    #if np.unique(obs['obi']) > 1:
    #    raise ValueError("Cannot determine roll for multi-obi obsid")
    # if there is more than one and they are close, just return the mean
    return np.mean(obs['roll'])


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


def get_fov(obsid, dbh):
    """
    For a given obsid, retrieve the pieces necessary to build an obsvis FOV.

    :param obsid: requested obsid
    :param aca_dbh: handle for aca database
    :returns: field-of-view
    :rtype: dict
    """

    target = dbh.fetchone("select * from target where obsid = %d" % obsid)

    if not target:
        logger.warn("No entry for %d in target table" % obsid)
        return None

    roll = archive_roll(obsid, dbh)
    if roll is None and target['lts_lt_plan']:
        target_date = DateTime("%04d:%03d" % (target['lts_lt_plan'].year,
                                              target['lts_lt_plan'].dayofyear))
        roll = Ska.Sun.nominal_roll(target['ra'], target['dec'], target_date)
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
               showtarget=1,
               showopticalaxis=1,
               showaimpoint=1,
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

    sim = dbh.fetchone("select * from sim where obsid = %d" % obsid)
    if sim:
        arcmin_mm = 2.925
        # just set one or the other (simzmm or simz) or obsvis complains
        # fov['offsetsimzmm'] = sim['trans_offset']
        fov['offsetsimz'] = sim['trans_offset'] / arcmin_mm

    if (target['acisid']):
        fov['unselectedchips'] = 0
        fov['subarrays'] = None
        acis = dbh.fetchone("""select * from acisparam a, target t
                               where a.acisid=t.acisid
                               and a.acisid=%d""" % target['acisid'])
        i_ccd_map = {'ccdi%d_on' % i: 'I%d' % i for i in range(0, 4)}
        s_ccd_map = {'ccds%d_on' % i: 'S%d' % i for i in range(0, 6)}
        ccd_map = dict(i_ccd_map.items() + s_ccd_map.items())
        ccd_list = []
        for x in ccd_map:
            if acis[x] == 'Y':
                ccd_list.append('%s' % ccd_map[x])
            if re.match('O.+', acis[x]):
                ccd_list.append('(%s)' % ccd_map[x])
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
        hrc = dbh.fetchone("""select * from hrcparam h, target t
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
            'showtarget',
            'showopticalaxis',
            'showaimpoint',
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

    dbh = Ska.DBI.DBI(dbi='sybase', server='sqlsao',
                      user=username, passwd=passwd,
                      database='axafocat')

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
            sobs = get_obsids_seqnum(reqid, dbh)
            if sobs is not None:
                logger.info("\tobsids: %s" % ' '.join([str(x) for x in sobs]))
                obsids.extend(sobs)
            else:
                logger.info("\tNo obsids found for sequence")
        if len(str(reqid)) == 8:
            logger.info("building fovs for obsids for proposal %d" % reqid)
            sobs = get_obsids_propnum(reqid, dbh)
            if sobs is not None:
                logger.info("\tobsids: %s" % ' '.join([str(x) for x in sobs]))
                obsids.extend(sobs)
            else:
                logger.info("\tNo obsids found for proposal")

    for obsid in obsids:
        fov = get_fov(obsid, dbh)

        if not fov:
            continue

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        outfile = os.path.join(outdir, 'obs%05d.fov' % obsid)
        fovfile = open(outfile, 'w')
        fov_lines = fov_text(fov)
        fovfile.writelines('\n'.join(fov_lines))
        fovfile.write("\n")
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
