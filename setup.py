from distutils.core import setup, Extension

mym_module = Extension('mym', 
                    include_dirs=['/usr/local/include', 
                                  './src/inc'
                                ],
                    #libraries=['tcl83'],
                    library_dirs=['/usr/local/lib'],
                    sources=['./mymmodule.c', 
                             './src/matrix.c',
                             './src/direct.c',
                             './src/iterative.c',
                             './src/eigenvalue.c',
                             './src/norm.c'
                            ],
                    #extra_compile_args="-o mym.so"
                            
                    )

setup(name="mym", 
    version="1.0", 
    description="this is a matrix library",
    author="wgtdkp",
    author_email="wgtdkp@sjtu.edu.cn",
    url="https://github.com/wgtdkp/mym",
    ext_modules=[mym_module])