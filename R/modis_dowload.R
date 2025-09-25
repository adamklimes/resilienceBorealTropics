## modis_download
library(MODIStsp)
# MODIStsp_get_prodnames()
# MODIStsp_get_prodlayers("MCD12Q1")

# credentials
# user <- [add username]
# pass <- [add password]
# SA (set start/end date)
MODIStsp(gui             = FALSE,
         out_folder      = paste0(getwd(), "/data/modis/SA"),
         out_folder_mod  = paste0(getwd(), "/data/modis/SA/HDFs"),
         selprod         = "Vegetation_Indexes_16Days_1Km (M*D13A2)",
         prod_version    = "061",
         bandsel         = "EVI",
         quality_bandsel = "QA_qual",
         sensor          = "Terra",
         user            = user,
         password        = pass,
         download_range  = "Seasonal",
         start_date      = "2017.12.01", # 2000
         end_date        = "2023.02.28", # 2006
         verbose         = TRUE,
         start_x         = 11,
         end_x           = 13,
         start_y         = 9,
         end_y           = 11,
         out_projsel     = "User Defined",
         output_proj     = "+proj=longlat +datum=WGS84 +no_defs",
         delete_hdf      = FALSE,
         parallel        = TRUE
         )

# RU
MODIStsp(gui             = FALSE,
         out_folder      = paste0(getwd(), "/data/modis/RU"),
         out_folder_mod  = paste0(getwd(), "/data/modis/RU/HDFs"),
         selprod         = "Vegetation_Indexes_16Days_1Km (M*D13A2)",
         prod_version    = "061",
         bandsel         = "EVI",
         quality_bandsel = "QA_qual",
         sensor          = "Terra",
         user            = user,
         password        = pass,
         download_range  = "Seasonal",
         start_date      = "2017.06.01", # 2000
         end_date        = "2022.08.31", # 2005
         verbose         = TRUE,
         start_x         = 19,
         end_x           = 23,
         start_y         = 1,
         end_y           = 4,
         out_projsel     = "User Defined",
         output_proj     = "+proj=longlat +datum=WGS84 +no_defs",
         delete_hdf      = FALSE,
         parallel        = TRUE
         )

# add rest of the year
# SA (set start/end date)
MODIStsp(gui             = FALSE,
         out_folder      = paste0(getwd(), "/data/modis/SA2"),
         out_folder_mod  = paste0(getwd(), "/data/modis/SA2/HDFs"),
         selprod         = "Vegetation_Indexes_16Days_1Km (M*D13A2)",
         prod_version    = "061",
         bandsel         = "EVI",
         quality_bandsel = "QA_qual",
         sensor          = "Terra",
         user            = user,
         password        = pass,
         download_range  = "Seasonal",
         start_date      = "2018.03.01", # 2001
         end_date        = "2022.11.30", # 2005
         verbose         = TRUE,
         start_x         = 11,
         end_x           = 13,
         start_y         = 9,
         end_y           = 11,
         out_projsel     = "User Defined",
         output_proj     = "+proj=longlat +datum=WGS84 +no_defs",
         delete_hdf      = FALSE,
         parallel        = TRUE
)

# RU
MODIStsp(gui             = FALSE,
         out_folder      = paste0(getwd(), "/data/modis/RU2"),
         out_folder_mod  = paste0(getwd(), "/data/modis/RU2/HDFs"),
         selprod         = "Vegetation_Indexes_16Days_1Km (M*D13A2)",
         prod_version    = "061",
         bandsel         = "EVI",
         quality_bandsel = "QA_qual",
         sensor          = "Terra",
         user            = user,
         password        = pass,
         download_range  = "Seasonal",
         start_date      = "2017.09.01", # 2000
         end_date        = "2022.05.31", # 2005
         verbose         = TRUE,
         start_x         = 19,
         end_x           = 23,
         start_y         = 1,
         end_y           = 4,
         out_projsel     = "User Defined",
         output_proj     = "+proj=longlat +datum=WGS84 +no_defs",
         delete_hdf      = FALSE,
         parallel        = TRUE
)


## land cover
# RU
MODIStsp(gui             = FALSE,
         out_folder      = paste0(getwd(), "/data/modis/RU_LC2"),
         out_folder_mod  = paste0(getwd(), "/data/modis/RU_LC2/HDFs"),
         selprod         = "LandCover_Type_Yearly_500m (MCD12Q1)",
         prod_version    = "006",
         bandsel         = "LC2",
         sensor          = "Terra",
         user            = user,
         password        = pass,
         download_range  = "Full",
         start_date      = "2018.01.01", # 2001
         end_date        = "2021.12.31", # 2003
         verbose         = TRUE,
         start_x         = 19,
         end_x           = 23,
         start_y         = 1,
         end_y           = 4,
         out_projsel     = "User Defined",
         output_proj     = "+proj=longlat +datum=WGS84 +no_defs",
         delete_hdf      = FALSE,
         parallel        = TRUE
)

# SA
MODIStsp(gui             = FALSE,
         out_folder      = paste0(getwd(), "/data/modis/SA_LC2"),
         out_folder_mod  = paste0(getwd(), "/data/modis/SA_LC2/HDFs"),
         selprod         = "LandCover_Type_Yearly_500m (MCD12Q1)",
         prod_version    = "006",
         bandsel         = "LC2",
         sensor          = "Terra",
         user            = user,
         password        = pass,
         download_range  = "Full",
         start_date      = "2018.01.01", # 2001
         end_date        = "2021.12.31", # 2003
         verbose         = TRUE,
         start_x         = 11,
         end_x           = 13,
         start_y         = 9,
         end_y           = 11,
         out_projsel     = "User Defined",
         output_proj     = "+proj=longlat +datum=WGS84 +no_defs",
         delete_hdf      = FALSE,
         parallel        = TRUE
)
#_
