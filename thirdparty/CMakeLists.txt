
if(CGOGN_USE_QT)
	add_subdirectory(libQGLViewer)
	set_target_properties(QOGLViewer PROPERTIES FOLDER thirdparty)
endif(CGOGN_USE_QT)

add_subdirectory(ply)
#add_subdirectory(OffBinConverter)
add_subdirectory(libMeshb)
add_subdirectory(termcolor)
