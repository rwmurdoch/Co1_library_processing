#!/bin/bash

find reads -type f -print0 | xargs -0 mv -t reads
