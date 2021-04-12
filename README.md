## 运行示例

1. 下载blast+，添加到环境变量

   ```bash
   cd ~/biosoft
   wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz
   tar -xzvf ncbi-blast-2.11.0+-x64-linux.tar.gz
   mv ncbi-blast-2.11.0+-x64-linux blast-2.11.0+
   export PATH=$PATH:~/biosoft/blast-2.11.0+/bin/
   
   makeblastdb -h
   ```

2. 将该项目克隆到本地

   ```
   git clone XXX
   ```

3. 运行代码执行测试

   ```bash
   # bash run.sh ErmC的序列 数据库文件
   # ErmC.fa 是用于在数据库中查找改序列的序列文件
   # 数据库文件，如果有多个数据文件
   # 例如：1.fa 2.fa 3.fa
   # 需要将它们合并成一个文件
   # cat 1.fa 2.fa 3.fa >database.fa
   bash run.sh ./test/ErmC.fa database.fa
   ```

4. 查看文件

   在Result文件夹下有两个文件

   + .upstream.nr.fa：非冗余的上游序列
   + .upstream.discard.fa：可能有问题的上游序列（与其他序列可能不同）

可能.upstream.discard.fa包含有用的序列，刚进行序列对齐建议将两个文件合并到一起。

```
cat *.upstream.nr.fa *.upstream.discard.fa > all.fa
```

认为.upstream.discard.fa包含的序列有误，就只采用.upstream.nr.fa文件。

## 其他

在`conf.sh`中有几个设置，自己查看里面的描述。

## 进一步

Mega比对或者CLUSTALW比对之后进行调整，然后可视化。
