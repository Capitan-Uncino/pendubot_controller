#!/bin/bash

nix-shell --run "nvidia-offload octave --quiet --eval 'wrapper("sim_ode45.m")'"
