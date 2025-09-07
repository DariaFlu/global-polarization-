file(REMOVE_RECURSE
  "libgp_dict.pdb"
  "libgp_dict.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/gp_dict.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
