from cook import create_task

# Create requirements task
create_task("requirements", action="pip-compile -v", targets=["requirements.txt"],
            dependencies=["requirements.in"])

# Create pip-sync task
create_task("pip-sync", action="pip-sync", dependencies=["requirements.txt"])