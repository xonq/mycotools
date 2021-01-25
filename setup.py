import setuptools

with open( "README.md", "r" ) as fh:
    long_description = fh.read()

setuptools.setup(
    name = "mycotools",
    version = "0.4.1b4",
    author = "xonq",
    author_email = "konkelzach@protonmail.com",
    description = "A compilation of bioinformatic and computation biology inspired " + \
        "python tools",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://gitlab.com/xonq/mycotools/mycotools",
    packages = setuptools.find_packages( ),
    scripts = ['mycotools/abstractHmm.py', 'mycotools/acc2fa.py', 'mycotools/acc2gff.py', 'mycotools/aggClus.py', 'mycotools/annotationStats.py', 'mycotools/assemblyStats.py', 'mycotools/curAnnotation.py', 'mycotools/db2hmmsearch.py', 'mycotools/fa2tree.py', 'mycotools/grabClusters.py', 'mycotools/jgiDwnld.py', 'mycotools/ncbiDwnld.py', 'mycotools/ome2name.py', 'mycotools/proteomeStats.py', 'mycotools/db2blast.py', 'mycotools/db/queryDB.py', 'mycotools/db/ncbi2db.py', 'mycotools/db/abstractDB.py', 'mycotools/db/predb2db.py', 'mycotools/db/updateDB.py', 'mycotools/db/jgi2db.py', 'mycotools/db/dbFiles.py', 'mycotools/db/gff2gff3.py'],
    install_requires = [ 'biopython', 'pandas', 'requests' ],
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires = '>=3.0,<4'
)
