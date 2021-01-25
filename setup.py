import setuptools

with open( "README.md", "r" ) as fh:
    long_description = fh.read()

setuptools.setup(
    name = "mycotools",
    version = "0.5.0b3",
    author = "xonq",
    author_email = "konkelzach@protonmail.com",
    description = "A compilation of bioinformatic and computation biology inspired " + \
        "python tools",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://gitlab.com/xonq/mycotools/mycotools",
    packages = setuptools.find_packages( ),
    scripts = ['mycotools/bin/abstractHmm.py', 'mycotools/bin/acc2fa.py', 'mycotools/bin/acc2gff.py', 'mycotools/bin/aggClus.py', 'mycotools/bin/annotationStats.py', 'mycotools/bin/assemblyStats.py', 'mycotools/bin/db2blast.py', 'mycotools/bin/db2hmmsearch.py', 'mycotools/bin/curAnnotation.py', 'mycotools/bin/jgiDwnld.py', 'mycotools/bin/ncbiDwnld.py', 'mycotools/bin/ome2name.py', 'mycotools/bin/fa2tree.py', 'mycotools/bin/gff2protein.py', 'mycotools/bin/grabClusters.py', 'mycotools/bin/proteomeStats.py', 'mycotools/bin/dbFiles.py', 'mycotools/utils/gff2gff3.py', 'mycotools/utils/jgi2db.py', 'mycotools/utils/ncbi2db.py', 'mycotools/utils/predb2db.py', 'mycotools/utils/abstractDB.py', 'mycotools/utils/queryDB.py', 'mycotools/utils/updateDB.py'],
    install_requires = [ 'biopython', 'pandas', 'requests' ],
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires = '>=3.0,<4'
)
