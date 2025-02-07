{
  "script_http"           : "/or.water.usgs.gov/proj/columbia",
  "script_cgi_bin"        : "/cgi-bin/columbia",
  "rasters"               : [ "tiffs/obtop.tif", "tiffs/smtop.tif", "tiffs/wntop.tif", "tiffs/grtop.tif", "tiffs/pmtop.tif"],
  "color_file"            : "framework_color_map.txt",
  "latlong_projection"    : "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", 
  "framework_projection"  : "+title=lcc +proj=lcc +lat_0=45.33333333333334 +lon_0=-120.5 +lat_1=47.33333333333334 +lat_2=45.83333333333334 +x_0=500000.0001016 +y_0=0 +ellps=GRS80 +datum=NAD83 +to_meter=0.3048006096012192 +no_defs",
  "northwest_x"   :      1210631.621,
  "northwest_y"   :      1189639.234,
  "northeast_x"   :      2951131.621,
  "northeast_y"   :      1189639.234,
  "southeast_x"   :      2951131.621,
  "southeast_y"   :      -405860.766,
  "southwest_x"   :      1210631.621,
  "southwest_y"   :      -405860.766,
  "noDataValue"   :  -32768,
  "xy_multiplier" :           0.3048,
  "xy_units"      :      "meters",
  "z_multiplier"  :            1.0,
  "z_units"       :      "feet",
  "graph_x_axis"  : "Distance, in feet",
  "graph_y_axis"  : "Elevation above NAVD 1988, in feet",
  "aboutFiles"    : {
                     "framework_welcome_txt"              : "framework_welcome.txt",
                     "framework_general_instructions_txt" : "framework_general_instructons.txt",
                     "framework_explanation_txt"          : "framework_explanation.txt",
                     "framework_contacts_txt"             : "framework_contacts.txt"
  },
  "wellLogFiles"    : {
                     "welllog_welcome_txt"                : "welllog_welcome.txt",
                     "framework_contacts_txt"             : "framework_contacts.txt"
  },
  "xsecFiles"    : {
                     "xsec_welcome_txt"                   : "xsec_welcome.txt",
                     "framework_contacts_txt"             : "framework_contacts.txt"
                    }
}
