
In order to run the batch system script 

1. list the datasets you need to use in .txt file with dataset names
1. change the workarea,exearea and X509 on to of the gefilelist.py 
1. run this command 

```bash
python getfilelist.py datesets.txt --out /path/to/out/dir
```
1. this will produce a .sh script called submit_all_to_batch_HTC.sh 

```bash
./submit_all_to_batch_HTC.sh 
```
it will do the job for you 


