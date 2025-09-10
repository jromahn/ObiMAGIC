#!/usr/bin/env perl

use Cwd 'abs_path';
use File::Basename;
use  File::Find;
use IO::Zlib;  #requires Compress::Zlib, of course
use IPC::Cmd qw[can_run run];
use List::MoreUtils qw(uniq);
