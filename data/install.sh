#!/bin/bash
echo "Downloading Remeff databases from Figshare..."
FIGSHARE_URL="https://figshare.com/ndownloader/files/66183590" #from https://doi.org/10.6084/m9.figshare.32859416
wget -O SUSSdata.tar.gz "$FIGSHARE_URL"
echo "Extracting databases..."
tar -xzf SUSSdata.tar.gz
rm SUSSdata.tar.gz
echo "Installation complete!"
