all:    main.o load_data.o maxwell.o
		nvcc -o all main.o load_data.o maxwell.o
 
main.o: main.cu
		nvcc -c main.cu -g

load_data.o: load_data.cu
		nvcc -c load_data.cu -g

maxwell.o: maxwell.cu
		nvcc -c maxwell.cu -g

clean:
		rm *.o all 
