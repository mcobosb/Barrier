#!/bin/bash

# Build and deploy script for Barrier project
# Cleans, builds, and copies executable to test directories

echo "=========================================="
echo "Barrier - Build and Deploy"
echo "=========================================="

# Clean previous build
echo "Cleaning previous build..."
make clean

# Build project
echo "Building project..."
cmake -DCMAKE_BUILD_TYPE=ReleaseFast -B .
make -C .

# Check if build was successful
if [ $? -eq 0 ]; then
    echo "Build successful!"
    
    # Copy to test directories
    echo "Deploying to test directories..."
    cp barrier tests/1022/barrier
    cp barrier tests/reducedGRE/barrier
    cp barrier tests/reducedGRE_z0/barrier
    cp barrier tests/qGRE/barrier
    cp barrier tests/x2/barrier
    
    echo "=========================================="
    echo "Deployment complete!"
    echo "  - tests/1022/barrier"
    echo "  - tests/reducedGRE/barrier"
    echo "  - tests/reducedGRE_z0/barrier"
    echo "  - tests/qGRE/barrier"
    echo "=========================================="
else
    echo "=========================================="
    echo "Build failed. Deployment skipped."
    echo "=========================================="
    exit 1
fi
