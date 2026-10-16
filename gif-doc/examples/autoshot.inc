
// ===== automatic frames (inserted by build_see_auto.sh, for the GIFs of the examples)
// Renders the conf-files from the one given on the command line, writes
// frame0000.tga, frame0001.tga, ... in the current folder, then exits.
// Environment: AUTOSHOT_STEP (conf-files between two frames, default 1),
// AUTOSHOT_LAST (last conf-file), AUTOSHOT_CMIN and AUTOSHOT_CMAX (fixed
// colour range for the whole animation).
static void autoShot(int) {
  int step = getenv("AUTOSHOT_STEP") ? atoi(getenv("AUTOSHOT_STEP")) : 1;
  int last = getenv("AUTOSHOT_LAST") ? atoi(getenv("AUTOSHOT_LAST")) : 1000000;
  glutReshapeWindow(width, height);
  if (getenv("AUTOSHOT_CMIN") && getenv("AUTOSHOT_CMAX")) {
    params["colorRangeMin"] = atof(getenv("AUTOSHOT_CMIN"));
    params["colorRangeMax"] = atof(getenv("AUTOSHOT_CMAX"));
    params["rescaleColorRange"] = 0;
    resetColors(params["colorMode"].get<int>(), 0);
  }
  int frame = 0;
  for (int k = 0; k < 3; k++) { textZone.reset(); display(); }  // let the window settle
  do {
    textZone.reset();
    display();  // twice: the screenshot may read the buffer of the previous frame
    display();
    char name[256];
    snprintf(name, 256, "frame%04d.tga", frame++);
    screenshot(name);
  } while (confNum + step <= last && tryToReadConf(confNum + step));
  exit(0);
}

