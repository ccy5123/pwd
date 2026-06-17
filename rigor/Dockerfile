# sci-adk-partition: the Docker experiment image for the partition-model capability.
#
# The stock sci-adk-python-base image carries numpy/scipy/networkx but NOT
# pandas/scikit-learn, which the partition CV experiment needs. This image adds them,
# pinned so a containerized run reproduces the committed CV numbers bit-for-bit.
#
# Build (on a Docker-enabled host, e.g. your local machine):
#     docker build -t sci-adk-partition rigor/
# Then `python rigor/run_model.py` auto-detects the daemon + image and runs the CV
# inside the container, stamping the image id into Evidence provenance.
FROM python:3.11-slim

WORKDIR /workspace

# numpy/pandas/scikit-learn install from prebuilt manylinux wheels on linux/amd64,
# so no apt build toolchain is needed (keeps the image small and the build hermetic).
COPY requirements-docker.txt /tmp/requirements-docker.txt
RUN pip install --no-cache-dir -r /tmp/requirements-docker.txt

CMD ["python3"]
