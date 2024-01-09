c This program creates a GrADS .gdat file that covers your
c   simulation domain.  This is a rectangular .gdat version of
c   the MicroMet mm.dat input file, where the rectancle covers
c   the largest extent of the mm points.  The .gdat data will
c   still be on the lat-lon (re)analysis grid, not the SnowModel
c   grid (it just, kind of, covers the SnowModel grid, plus the
c   (re)analysis grid cells that are just outside the SnowModel
c   domain boundary.  This file does two things: it allows you
c   a nice way to look at the met forcing that is going into the
c   SnowModel simulation, and it is used to correct (re)analysis
c   drizzle (see the 5_fix_drizzle_optional/ directory).
