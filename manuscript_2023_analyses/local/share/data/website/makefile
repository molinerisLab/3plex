#to use Jekyll work on ncrawler

experiment_table.csv:
	https://docs.google.com/spreadsheets/d/1-srAK7JpCCSuYx6T8U-aYQW9xuqt3pJ2BPBMTIOZNNk/edit#gid=0

experiment_table.html: experiment_table.csv
	@echo consider to conda acrivate py3env
	python -c '\
		import pandas as pd; \
		data=pd.read_csv("$<");\
		data_html=data.to_html();\
		out_file=open("$@","w");\
		out_file.write(data_html)\
	'
