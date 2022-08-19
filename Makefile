hbv: main.c arrays.c hbv_model.c hbv_report.c hbv_performance.c readcsv.c 
	gcc -o hbv main.c arrays.c hbv_model.c hbv_report.c hbv_performance.c readcsv.c -lm -fopenmp -Wall

#hbvgpu: main.c arrays.c hbv_model.c hbv_report.c hbv_performance.c readcsv.c 
#	gcc -o hbv main.c arrays.c hbv_model.c hbv_report.c hbv_performance.c readcsv.c -lm -fopenmp -foffload=nvptx-none -Wall

clean:
	rm -f hbv
