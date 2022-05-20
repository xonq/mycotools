import setuptools

with open( "README.md", "r" ) as fh:
    long_description = fh.read()

setuptools.setup(
    name = "mycotools",
    version = "0.11b11",
    author = "xonq",
    author_email = "konkelzach@protonmail.com",
    description = "A compilation of bioinformatic and computation biology inspired " + \
        "python tools",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://gitlab.com/xonq/mycotools/mycotools",
    package_dir={"": "mycotools"},
    packages = setuptools.find_packages( where="mycotools" ),
    scripts = ['mycotools/mycotools/acc2fa.py', 'mycotools/mycotools/acc2gff.py', 'mycotools/mycotools/aggClus.py', 'mycotools/mycotools/annotationStats.py', 'mycotools/mycotools/assemblyStats.py', 'mycotools/mycotools/extractDB.py', 'mycotools/mycotools/curAnnotation.py', 'mycotools/mycotools/db2hmm.py', 'mycotools/mycotools/dbFiles.py', 'mycotools/mycotools/fa2tree.py', 'mycotools/mycotools/gff2svg.py', 'mycotools/mycotools/jgiDwnld.py', 'mycotools/mycotools/ncbiDwnld.py', 'mycotools/mycotools/ome2name.py', 'mycotools/mycotools/acc2loci.py', 'mycotools/mycotools/predb2db.py', 'mycotools/mycotools/db2search.py', 'mycotools/mycotools/crap.py', 'mycotools/mycotools/fa2hmmer2fa.py', 'mycotools/mycotools/fa2weight.py', 'mycotools/mycotools/gff2seq.py', 'mycotools/mycotools/coords2fa.py', 'mycotools/mycotools/mycodb', 'mycotools/mycotools/utils/queryDB.py'],
    install_requires = [ 'biopython', 'pandas', 'requests', 'dna-features-viewer', 'scipy' ],
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires = '>=3.0,<4'
)
