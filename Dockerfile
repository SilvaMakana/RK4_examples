# Use an official GCC runtime as a parent image
FROM gcc:14-bookworm

# Install gnuplot
RUN apt-get update && \
    apt-get install -y gnuplot gnuplot-nox && \
    rm -rf /var/lib/apt/lists/*

# Set the working directory in the container
WORKDIR /workspace/RK4_Kepler

# Copy the local C source code to the container
COPY . .

