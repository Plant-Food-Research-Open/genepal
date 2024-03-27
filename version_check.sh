#!/usr/bin/env bash

config_version=$(sed -n "s/.*version.*= '\(.*\)'.*/\1/p" nextflow.config)

# Check CHANGELOG version

grep "## $config_version - " CHANGELOG.md >/dev/null \
    || (echo 'Failed to match CHANGELOG version'; exit 1)
