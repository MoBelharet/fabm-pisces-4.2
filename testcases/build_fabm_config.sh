#*************************************************************************************
# Author : Mokrane BELHARET
# Company: Mercator Ocean international
# date : October 2023
#------------------------------
# Description:
# This schell script helps to set-up and compile the fabmized biogeochemical models
#************************************************************************************


# Delete the build folder if present
rm -rf build
#************* PARAMETERS *******************************

# The project directory 
project_directory=`pwd`

# name of the folder in which cmake files will be saved after compilation
build_folder="build/fabm"

# path to fabm folder
path_to_fabm="/home/mbelharet/FABM_PROJECT/fabm"      #"/home/mbelharet/FABM_PROJECT/PROJECT_NEMO4.0_FABM_PISCES/fabm"

# The name of the physical model to be used
hostname="nemo"

#Biogeochemical model
institute="pisces"

#Other models
other_models=true
models=migration  # add all the other models separated by comma . ex: toto1,toto2,toto3
   
#Compiler
change_compiler=true
compiler="mpif90" #To add only if change_compiler is true . gfortran is the compiler by defauld

#Path where the compiled code will create fabm library and necessary include files (By default it is : ~/local/fabm/$hostname)
change_compiled_fabm_path=true
path_compiled_fabm="$project_directory/local/fabm/$hostname" #To add only if change_compiled_fabm_path is true

# Activate Dbug mode
activate_debug=false


#************* SETTING UP  ****************
# PLEASE DON'T MODIFY THE FOLLOWING PARTS
#****************************************** 
config="$path_to_fabm -DFABM_HOST=$hostname"


path_to_model="$project_directory/fabm-$institute" 
DFABM_model_BASE="-DFABM_${institute^^}_BASE=${project_directory}/fabm-${institute}"

if [ $other_models = true ]; then
  for model in ${models//,/ }; do
    institute="${institute};${model}"
    
    DFABM_model_BASE="$DFABM_model_BASE -DFABM_${model^^}_BASE=${project_directory}/fabm-${model}" 
  done
fi

config="$config -DFABM_INSTITUTES=$institute ${DFABM_model_BASE}"

if [ $change_compiler = true ]; then
	config="$config -DCMAKE_Fortran_COMPILER=$compiler"
fi

if [ $change_compiled_fabm_path = true ]; then
	config="$config -DCMAKE_INSTALL_PREFIX=$path_compiled_fabm"
fi

if [ $activate_debug = true ]; then
	config="$config -DCMAKE_BUILD_TYPE=Debug" 
fi


echo $config

#************ BUILDING ************************

mkdir -p $build_folder
cd $build_folder

cmake $config
make install

cd $project_directory


