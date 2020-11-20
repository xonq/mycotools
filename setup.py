import setuptools

with open( "README.md", "r" ) as fh:
	long_description = fh.read()

setuptools.setup(
	name = "kontools-xonq",
	version = "0.0.1",
	author = "xonq",
	author_email = "konkelzach@protonmail.com",
	description = "A compilation of bioinformatic and computation biology inspired " + \
		"python tools",
	long_description = long_description,
	long_description_content_type = "text/markdown",
	url = "https://gitlab.com/xonq/kontools"
	packages = setuptools.find_packages(),
	classifiers = [
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: Linux / Unix",
	],
	python_requires = '>=3.0',
)
