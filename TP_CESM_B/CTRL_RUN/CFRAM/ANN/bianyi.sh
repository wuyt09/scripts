#!/bin/bash
ifort baseline.f -o base.out &
ifort warm.f -o warm.out &
ifort drdt.f -o drdt.out &
ifort albedo.f -o albedo.out &
ifort co2.f -o co2.out &
ifort o3.f -o o3.out &
ifort solar.f -o solar.out &
ifort t.f -o t.out &
ifort wv.f -o wv.out &
ifort cloud.f -o cloud.out &
ifort cfram_ranc.f -o cfram.out &
