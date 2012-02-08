#!/usr/bin/env python

import re
import os
import getpass
import urllib
from BeautifulSoup import BeautifulSoup

import asciitable
import numpy as np
import Ska.DBI

from Chandra.Time import DateTime

import logging
logger = logging.getLogger()

last_aimpoint_shift = DateTime('2011:187')

def lts_table(lts_url='http://cxc.harvard.edu/target_lists/longsched.html'):
    soup = BeautifulSoup(urllib.urlopen(lts_url).read())
    page = "".join(soup.findAll(text=True)) 
    pagelines = page.splitlines()
    obslines = []
    for line in pagelines:
        # the lines that begin with a sequence number and obsid and end with
        # proposal number
        if re.match('^\d{6}\s+(P1)?\s+\d{1,5}.*\d{8}\s*$', line):
            obslines.append(line)

    col_starts=(0, 7, 10, 16, 38, 45, 53, 61, 67, 
                73, 80, 82, 84, 89, 110, 115, 
                119, 122, 125, 128, 131, 134, 
                138, 140, 165)
    col_ends=(6, 9, 14, 36, 43, 52, 60, 66, 72, 
              79, 81, 83, 88, 107, 112, 117, 
              120, 123, 126, 129, 132, 135, 
              139, 162, 173)
    names=('seqnbr', 'pool', 'obsid', 'target', 'ksec', 'ra', 'dec', 'roll', 'pitch',
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

def prelim_roll(obsid, lts=lts_table()):
    match = lts['obsid'] == obsid
    n_match = len(np.flatnonzero(match))
    if n_match == 0:
        return None
    if n_match > 1:
        raise ValueError
    return lts[match]['roll'][0]

def archive_roll(obsid, dbh):
    obs = dbh.fetchall("""select obsid, obi, ra_targ, dec_targ, ra_pnt, dec_pnt, roll_pnt
from observations where obsid = %d""" % obsid)
    if not len(obs):
        return None
    if np.unique(obs['obi']) > 1:
        raise ValueError("Cannot determine roll for multi-obi obsid")
    if np.std(obs['roll_pnt'] > 2):
        raise ValueError("Observed rolls for obsid have std > 2.00; failing")
    return np.mean(obs['roll_pnt'])

def get_obsids_seqnum(reqid, dbh):
    targs = dbh.fetchall("select obsid from target where seq_nbr = '%d'" % reqid)
    if len(targs):
        return targs['obsid']

def get_obsids_propnum(reqid, dbh):
    targs = dbh.fetchall("""select obsid from target t 
join prop_info p on t.ocat_propid = p.ocat_propid 
where prop_num = '%d'""" % reqid);
    if len(targs):
        return targs['obsid']

def default_offsets(det):
    defaults = {'HRC-I': (0.0,0.0),
                'HRC-S': (0.0,0.0),
                'ACIS-I': (-0.2, -0.25),
                'ACIS-S': (0.15, -0.25)}
    return defaults[det]


def main(reqids=[], outdir='.'):

    console = logging.StreamHandler()
    logger.setLevel(logging.INFO)
    logger.addHandler(console)

    arcmin_mm = 2.925
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


    username = getpass.getuser()
    dbpassfile = os.path.join(os.environ['HOME'], '.arc5gl_pwd')

    if os.path.exists(dbpassfile):
        passwd = open(dbpassfile).read().strip()
    else:
        passwd = getpass.getpass()

    sqlsao = Ska.DBI.DBI(dbi='sybase', server='sqlocc', 
                      user=username, passwd=passwd,
                      database='axafocat')
    sqlaca = Ska.DBI.DBI(dbi='sybase', server='sybase',
                         user='aca_read', database='aca')

    obsids = []
    for reqid in reqids:
        if len(str(reqid)) <= 5:
            logger.info("building fovs for obsid %d" % reqid)
            obsids.extend([reqid])
        if len(str(reqid)) == 6:
            logger.info("building fovs for obsids for seqnum %d" % reqid)
            sobs = get_obsids_seqnum(reqid, sqlsao)
            logger.info("\tobsids: %s" % ' '.join([str(x) for x in sobs]))
            obsids.extend(sobs)
        if len(str(reqid)) == 8:
            logger.info("building fovs for obsids for proposal %d" % reqid)
            sobs = get_obsids_propnum(reqid, sqlsao)
            logger.info("\tobsids: %s" % ' '.join([str(x) for x in sobs]))
            obsids.extend(sobs)

    for obsid in obsids:

        target = sqlsao.fetchone("select * from target where obsid = %d" % obsid)        

        roll = archive_roll(obsid, sqlaca)
        if roll is None:
            lts = lts_table()
            roll = prelim_roll(obsid, lts)
        if roll is None:
            roll = 0

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

        if (fov['offsety'] is None) or (fov['offsetz'] is None):
            def_offsets = default_offsets(target['instrument'])
            if fov['offsety'] is None:
                fov['offsety'] = def_offsets[0]
            if fov['offsetz'] is None:
                fov['offsetz'] = def_offsets[1]

        sim = sqlsao.fetchone("select * from sim where obsid = %d" % obsid)
        if sim:
        #    fov['offsetsimzmm'] = sim['trans_offset']
            fov['offsetsimz'] = sim['trans_offset'] / arcmin_mm

        if (target['acisid']):
            fov['unselectedchips'] = 0
            fov['subarrays'] = None
            acis = sqlsao.fetchone("""select * from acisparam a, target t
                                   where a.acisid=t.acisid
                                   and a.acisid=%d""" % target['acisid'])
            i_ccd_map = {'ccdi%d_on' % i : 'I%d' % i for i in range(0,4)}
            s_ccd_map ={'ccds%d_on' % i : 'S%d' % i for i in range(0,6)}
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
            hrc = sqlsao.fetchone("""select * from hrcparam h, target t
                                     where h.hrcid=t.hrcid
                                     and h.hrcid=%d""" % target['hrcid'])
            if target['instrument'] == 'HRC-S':
                fov['unselectedchips'] = 0
                fov['timingmode'] = 0
                if hrc['timing_mode'] == 'Y':
                    fov['timingmode'] = 1
        

        if target['lts_lt_plan'] is not None:
            ot = target['lts_lt_plan']
            ot_date = DateTime('%04d:%03d' % (ot.year, ot.dayofyear))

        if target['soe_st_sched_date'] is not None:
            ot = target['soe_st_sched_date']
            ot_date = DateTime('%04d:%03d' % (ot.year, ot.dayofyear))

        if ot_date and ot_date.secs < last_aimpoint_shift.secs:
            logger.warn("inaccurate fov: obsid %d predates last aimpoint shift." % obsid)
                              
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outfile = os.path.join(outdir, 'obs%05d.fov' % obsid)
        fovfile = open(outfile, 'w')
        fovfile.write("fovid: %(fovid)s\n" % fov)
        for col in cols:
            if col in fov:
                fovfile.write("   %s: %s\n" % (col, fov[col]))
        fovfile.close()
        logger.info("obsid %d fov written to: %s" % (obsid, outfile))


def get_args():
    from argparse import ArgumentParser
    parser = ArgumentParser(description="obsvis fov generator from ocat")
    parser.add_argument('reqids', 
                        type=int,
                        nargs='+',
                        help="id or ids to fetch")
    parser.add_argument("--outdir", 
                        default=".",
                        help="output directory for fov files")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    my_args = get_args()
    main(**my_args.__dict__)



