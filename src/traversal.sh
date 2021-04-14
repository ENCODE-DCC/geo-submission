#!/bin/bash

while read name
do
    ascp -i /home/ubuntu/private.openssh -T --policy=fair -l600m -k1 "$name" asp-dcc-stanford@upload.ncbi.nlm.nih.gov:incoming
done
