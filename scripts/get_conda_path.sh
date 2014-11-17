#!/bin/bash
conda info --system | grep sys.prefix | cut -d: -f2 | sed 's/^ *//'

