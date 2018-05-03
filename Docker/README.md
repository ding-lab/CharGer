# Dockerfile for CharGer

Comments:

Currently CharGer supports the Python 2.x series. This Dockerfile does not build a local version of the Variant Effect Predictor (VEP) tool. Accordingly, CharGer will attempt to use Ensembl's REST API to supply VEP annotation if that is not already provided in the input variant call files.


How to build:

	docker build -t charger .
