void setup_xsec_analyzer() {

  int return_code = 0;

  // Load the xsec_analyzer shared library
  return_code = gSystem->Load("libXSecAnalyzer");
  if (return_code == 0) std::cout << "\nSuccessfully loaded XSecAnalyzer"
    " dynamic library.\n";
  else std::cout << "\nError loading XSecAnalyzer dynamic library.\n"
    << "Please add the directory that contains this library to your\n"
    << "LD_LIBRARY_PATH environment variable (or equivalent on\n"
    << "non-Linux systems) and try again.\n";
}
