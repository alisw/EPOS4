//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include <sys/times.h>
#include <sys/timeb.h>
#include <stdio.h>
#include <malloc.h>
#include <unistd.h>

/* usage from Fortran:  call timer(iutime)  */

void timer_(double *etime){
struct timeb tr;
struct tms tu;
ftime (&tr);
times (&tu);
etime[0]=(double) tu.tms_utime;
etime[1]=(double) tu.tms_stime;
etime[2]=(double) tr.time;
etime[3]=(double) tr.millitm;
etime[4]=(double) tr.timezone;
//printf("timer %f %f \n", etime[2],etime[3]);
}

void checkmemory_(int *memory)
{
        char buf[30];
        snprintf(buf, 30, "/proc/%u/statm", (unsigned)getpid());
        FILE* pf = fopen(buf, "r");
        if (pf) {
            unsigned size; //       total program size in pages (page=4kB)
            unsigned resident;//   resident set size
            unsigned share;//      shared pages
            unsigned text;//       text (code)
            unsigned lib;//        library
            unsigned data;//       data/stack
            unsigned dt;//         dirty pages (unused in Linux 2.6)
            fscanf(pf, " %u %u %u %u %u %u", &size, &resident, &share, &text, &lib, &data);
            /*printf(" memory = %u ",size);*/
            *memory=size*4/1000;   //MB
        }
        fclose(pf);
}        
