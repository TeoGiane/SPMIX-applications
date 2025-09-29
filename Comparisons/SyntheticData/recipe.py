from cook import create_task # , Task
# from cook.contexts import create_group
# import itertools as it
# import numpy as np
# import os
# from pathlib import Path
# from typing import Literal

# Create requirements task
create_task("requirements", action="pip-compile -v", targets=["requirements.txt"],
            dependencies=["requirements.in"])

# Create pip-sync task
create_task("pip-sync", action="pip-sync", dependencies=["requirements.txt"])

# Create generate-data task
action = ["Rscript","--vanilla","src/generate_data.R"]
deps = ["src/generate_data.R"]
create_task(name="generate_data", dependencies=deps, action=action)
