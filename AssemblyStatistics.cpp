#include "iostream"
#include "stdio.h"
#include "string"
#include "string.h"
#include "stdint.h"
#include "vector"
#include "cstdlib"
#include <map>
#include "memory"
#include <algorithm>
#include "fstream"
#include "sstream"
#include "list"
#include "stdlib.h"
#include "time.h"
using namespace std;

bool get_a_fasta_read(ifstream & fasta_in, string &tag, string &str, string & n_tag)
{

	ifstream tmp_ifstream;
	string temp;
	if (!getline(fasta_in, temp))
	{
		return 0;
	}
	if (temp[temp.size() - 1] == '\n' || temp[temp.size() - 1] == '\r')
	{
		temp.resize(temp.size() - 1);
	}

	str.clear();
	if (temp[0] == '>')
	{
		tag = temp;
	}
	else
	{
		tag = n_tag;
		str = temp;
	}


	while (getline(fasta_in, temp))
	{

		if (temp[temp.size() - 1] == '\n' || temp[temp.size() - 1] == '\r')
		{
			temp.resize(temp.size() - 1);
		}

		if ((temp.size()>0 && (temp[0] == '>' || temp[0] == '\n' || temp[0] == '\r')))
		{
			n_tag = temp;
			return 1;
		}
		else
		{
			str += temp;

		}

	}
	return 1;
}


bool get_a_fastq_read(ifstream & fastq_in, string &tag, string &seq, string & quality)
{

	ifstream tmp_ifstream;
	string temp;
	if (!getline(fastq_in, temp))
	{
		return 0;
	}
	seq.clear();
	if (temp[0] == '@')
	{
		tag = temp;
	}
	else
	{
		return 0;
	}
	getline(fastq_in, seq);//seq
	if (seq[seq.size() - 1] == '\n' || seq[seq.size() - 1] == '\r')
	{
		seq.resize(seq.size() - 1);
	}
	getline(fastq_in, temp);//'+'
	getline(fastq_in, quality);
	if (quality[quality.size() - 1] == '\n' || quality[quality.size() - 1] == '\r')
	{
		quality.resize(quality.size() - 1);
	}

	return 1;
}


int main(int argc, char* argv[])
{

	cout<<"Command:  ProgramFile contigs contigs_filename LenTh cut_off_length"<<endl; 
	string in_fname="LongContigs.txt";
	int LenTh=100;
	bool JustLengths=0;
	bool ContigMode=1;
	uint64_t GenomeSize=0;
	for(int i=1;i<argc;++i)
	{

		if(strcmp(argv[i],"info")==0)
		{
			i++;
			in_fname=(argv[i]);
			ContigMode=0;
			continue;
		}
	
		if(strcmp(argv[i],"contigs")==0)
		{
			i++;
			in_fname=(argv[i]);
			ContigMode=1;
			continue;
		}
		if(strcmp(argv[i],"GS")==0)
		{
			i++;
			sscanf(argv[i],"%Lu",&GenomeSize);
			//GenomeSize=atoi(argv[i]);
			
			continue;
		}
		if(strcmp(argv[i],"LenTh")==0)
		{
			i++;
			LenTh=atoi(argv[i]);
			continue;
		}
	}
	
	ifstream in_long_contigs(in_fname.c_str());
	string s,s1,s2,s3,s4,s5,s6,tag1,tag2,seq_s;
	vector<size_t> contig_len_vt;
	int contig_no;
	size_t cont_len;
	cout<<"Loading..."<<endl;
	getline(in_long_contigs, s1);
	bool fq_flag = 0;
	if (s1.size() == 0)
	{
		return -1;
	}
	
	if (s1[0] == '>' || s1[0] == '@')
	{
		if (s1[0] == '@')
		{
			fq_flag = 1;
		}
	}
	else
	{
		JustLengths = 1;
		ContigMode = 0;
	}
	in_long_contigs.close();
	in_long_contigs.open(in_fname.c_str());
	int nr = 0;
	string tag, tag_n, contig_str,quality;
	bool read_success = 1;
	if(JustLengths==0)
	{

		
		while(read_success)
		{

			if (fq_flag)
			{
				read_success = get_a_fastq_read(in_long_contigs, tag, contig_str, quality);
			}
			else
			{
				read_success = get_a_fasta_read(in_long_contigs, tag, contig_str, tag_n);
			}

			nr++;
			if (read_success)
			{
				contig_len_vt.push_back(contig_str.size());

			}
			
		}
	}
	else
	{
			
		while(in_long_contigs>>cont_len)
		{
			contig_len_vt.push_back(cont_len);
		}
		
	}

	
	
	cout<<"Sorting..."<<endl;
	sort(contig_len_vt.begin(),contig_len_vt.end());
	//vector<size_t> cvt2=contig_len_vt;
	//reverse(cvt2.begin(),cvt2.end());
	//ofstream o_temp("temp.txt");
	//for (int i=0;i<100;++i)
	//	o_temp<<cvt2[i]<<endl;
	uint64_t sum=0,sum1=0;
	for (size_t i=0;i<contig_len_vt.size();++i)
	{
		if(contig_len_vt[i]>LenTh)
		{
			sum+=contig_len_vt[i];
		}
	}

	string o_file=in_fname+"Stats.txt";
	ofstream o_stats(o_file.c_str());
	o_file=in_fname+"Stats_10k_5k_1k.txt";
	ofstream o_stats2(o_file.c_str());
	sum1=0;
	uint64_t sum2=0;
	double r=0.1;
	int n_stat=10;
	o_stats<<"Total len: "<<sum<<endl;
	o_stats<<"Longest contig: "<<contig_len_vt[contig_len_vt.size()-1]<<endl;
	if(GenomeSize>0)
	{
		sum=GenomeSize;
	}
	for (int i=contig_len_vt.size()-1;i>=1;--i)
	{
		sum1+=contig_len_vt[i];
	
		if(contig_len_vt[i-1]<10000&&contig_len_vt[i]>=10000)
			o_stats2<<"10 kbp: #"<<contig_len_vt.size()-i<<" ,add up to: "<<sum1<<endl;
		if(contig_len_vt[i-1]<5000&&contig_len_vt[i]>=5000)
			o_stats2<<"5 kbp: #"<<contig_len_vt.size()-i<<" ,add up to: "<<sum1<<endl;
		if(contig_len_vt[i-1]<1000&&contig_len_vt[i]>=1000)
			o_stats2<<"1 kbp: #"<<contig_len_vt.size()-i<<" ,add up to: "<<sum1<<endl;
		if(contig_len_vt[i-1]<100&&contig_len_vt[i]>=100)
		{
			o_stats2<<"100 bp: #"<<contig_len_vt.size()-i<<" ,add up to: "<<sum1<<endl;
			o_stats2<<"Avg contig sz: "<<(int)((double)sum1/double(contig_len_vt.size()-i)+0.5)<<endl;
		}


		while(sum1>double(sum)*r)
		{
			if(r<=0.91)
			{
			o_stats<<"N"<<n_stat<<" "<<contig_len_vt[i]<<endl;
			r+=0.1;
			n_stat+=10;
			}
			if(r>0.91)
			{break;}
		}
		
		
	}


}