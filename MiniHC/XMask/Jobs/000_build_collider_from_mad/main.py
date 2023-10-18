
"""This script is used to build the base collider with Xmask, configuring only the optics. Functions
in this script are called sequentially."""
# ==================================================================================================
# --- Imports
# ==================================================================================================
from cpymad.madx import Madx
import os
import xmask as xm
import xmask.lhc as xlhc
import shutil
import yaml
from pathlib import Path

import xtrack as xt
import xpart as xp




# ==================================================================================================
# --- Function to build collider from mad model
# ==================================================================================================
def build_collider_from_mad(madfile):
    # Importing LHC sequences:
    lines = {}
    for seq in ['lhcb1','lhcb2']:

        with open(f'madx_input_{seq}.log', 'w') as f:
            mad = Madx(stdout=f)
            mad.option(echo = True, warn = True)

            file_content = open(madfile).read()
            
            if seq == 'lhcb2':
                file_content = file_content.replace('lhcb1 : line= (clockwise);','lhcb2 : line= (anticlockwise);')
                file_content = file_content.replace('beam.bv.b1    = 1      ;',
                                                    'beam.bv.b2    = -1     ;')
                file_content = file_content.replace('b1','b2')
            # else:
            #     file_content = file_content.replace('px = xing.ip5   ,','px = -xing.ip5   ,')
            #     pass

            mad.input(file_content)
            mad.sequence[seq].use()


            lines[seq] = xt.Line.from_madx_sequence(    mad.sequence[seq],
                                                        deferred_expressions=True)
            
            mad_beam = mad.sequence[seq].beam
            lines[seq].particle_ref = xp.Particles(p0c = mad_beam.pc*1e9,q0 = mad_beam.charge, mass0 = mad_beam.mass*1e9)

    collider = xt.Multiline(lines=lines)


    return collider



def clean():
    # Remove all the temporaty files created in the process of building collider
    os.remove("madx_input_lhcb1.log")
    os.remove("madx_input_lhcb2.log")


# ==================================================================================================
# --- Main function for building distribution and collider
# ==================================================================================================
def build_collider(madfile = "../../../Backend/MAD/MiniHC.mad",
                   export_json = 'zfruits/collider_000.json',
                   sanity_checks = True):



    # Build collider from mad model
    collider = build_collider_from_mad(madfile=madfile)

    # Twiss to ensure eveyrthing is ok
    if sanity_checks:
        collider["lhcb1"].twiss(method="4d")
        collider["lhcb2"].twiss(method="4d")

    # Clean temporary files
    clean()


    collider.to_json(export_json)


# ==================================================================================================
# --- Script for execution
# ==================================================================================================

if __name__ == "__main__":
    build_collider()






