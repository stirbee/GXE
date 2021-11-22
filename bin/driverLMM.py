# -*- coding: utf-8 -*-
import multiprocessing,os,sys,datetime,gzip


def driver(n):

	os.system("python ../bin/GXEWAS.py ./cv%d/%s  ./cv%d/%s ./cv%d/%s  ./cv%d/%s  >& /dev/null"%(n, BED, n, PHE, n,ENV, n,OUT))
	#os.system("cp ../result.csv ./cv%d >& /dev/null"%(n))		


if __name__ == "__main__":
	print ('start time:', datetime.datetime.now().strftime("%Y/%m/%d %H:%S:%M"))
	print ('cpu number:',  multiprocessing.cpu_count())
	

	
	BED = sys.argv[1]
	PHE = sys.argv[2]
	ENV = sys.argv[3]
	OUT = sys.argv[4]
	N=int(sys.argv[5])  #########设置进程数
	pool = multiprocessing.Pool(N)
	
	
	for i in range(1,N+1):
		pool.apply_async(driver,(i,))   #维持执行的进程总数为processes，当一个进程执行完毕后会添加新的进程进去
	
	pool.close()
	pool.join()   #调用join之前，先调用close函数，否则会出错。执行完close后不会有新的进程加入到pool,join函数等待所有子进程结束
	
	print ("Sub-process(es) done.")
	print ('finished time:',datetime.datetime.now().strftime("%Y/%m/%d %H:%S:%M"))	
		
