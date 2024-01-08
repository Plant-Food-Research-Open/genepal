#!/usr/bin/env bash

nextflow \
    main.nf \
    -profile local,docker \
    -resume \
    -stub \
    --params-file conf/local_stub_params.json