#!/bin/bash

git pull --rebase origin master
git submodule update connectivity-diagnostics
git submodule update ggmClass
