#!/usr/bin/env bash

echo "Cleaning out intermediate files..." >&2
rm genome/*.fa.*
rm tmp/*
echo "Done." >&2
