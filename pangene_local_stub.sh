#!/usr/bin/env bash

nextflow \
    main.nf \
    -profile local,docker \
    -resume \
    -stub \
    --max_cpus=1 \
    --max_memory=1.GB