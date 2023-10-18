
from pathlib import Path
import importlib


main = importlib.import_module('MiniHC.XMask.Jobs.000_build_collider_from_mad.main')


if not Path('MiniHC/JSON').exists():
        Path('MiniHC/JSON').mkdir()

main.build_collider(madfile     = "MiniHC/MAD/MiniHC.mad",
                    export_json = "MiniHC/JSON/MiniHC.json",
                    sanity_checks = True)