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
	ofstream o_log("LongReadSelection_log.txt");
	cout << "Command:  ProgramFile sum total_length longest 0 o outfile f fa/fq_file f fa/fq_file " << endl;
	uint64_t sum = 0,total_sum=0;
	bool longest = 0;
	vector<int> read_len_vec;
	vector<string> files;
	vector<string> read_names;
	string filename_out="out.fasta";
	for (int i = 1; i<argc; ++i)
	{

		if (strcmp(argv[i], "f") == 0)
		{
			i++;
			files.push_back (argv[i]);
			
			continue;
		}
		if (strcmp(argv[i], "o") == 0)
		{
			i++;
			filename_out=(argv[i]);

			continue;
		}
		
		if (strcmp(argv[i], "sum") == 0)
		{
			i++;
			
			stringstream strValue;
			strValue << argv[i];
			strValue >> total_sum;

			continue;
		}
		if (strcmp(argv[i], "longest") == 0 || strcmp(argv[i], "Longest") == 0)
		{
			i++;
			longest = atoi(argv[i]);

			continue;
		}
	}
	sum = 0;
	
	ofstream o_selected_reads;
	if (!longest)
	{
		o_selected_reads.open(filename_out.c_str());
	}
	for (int i = 0; i < files.size(); ++i)
	{
		ifstream read_in(files[i].c_str());
		cout << "Loading file: " << files[i].c_str() << endl;
		string temp;
		getline(read_in, temp);
		bool fq_flag = 0;
		if (temp[0] == '@')
		{
			fq_flag = 1;
		}
		read_in.close();

		size_t n_reads = 0;
		string tag, n_tag,seq,qs;
		bool read_success = 1;
		read_in.open(files[i].c_str());
		while (read_success)
		{
			if (!longest && (sum > total_sum))
			{
				break;
			}
			if (fq_flag)
			{
				read_success = get_a_fastq_read(read_in,tag,seq,qs);
			}
			else
			{
				read_success = get_a_fasta_read(read_in, tag, seq,n_tag);
			}
			if (read_success)
			{
				n_reads++;
			}
			else
			{
				break;
			}
			if (!longest)
			{
				tag[0] = '>';
				o_selected_reads << tag << endl;
				o_selected_reads << seq << endl;
			}

			read_len_vec.push_back(seq.size());
			sum += seq.size();
			
		}
		cout << n_reads << " reads loaded." << endl;
		
		read_in.close();
	}
	cout << "Total bases: " << sum << endl;
	o_log << "Total bases: " << sum << endl;
	if (!longest)
	{
		return 0;
	}
	sum = 0;
	cout << "Sorting..." << endl;
	sort(read_len_vec.begin(), read_len_vec.end());

	int trim_len = 0;
	for (vector<int>::reverse_iterator rit = read_len_vec.rbegin(); rit != read_len_vec.rend(); ++rit)
	{
		sum += *rit;
		if (sum > total_sum)
		{
			trim_len = *rit;
			break;
		}
	}
	cout << "Trim length: " << trim_len << endl;
	o_log << "Trim length: " << trim_len << endl;
	size_t n_reads = 0;
	sum = 0;
	o_selected_reads.open(filename_out.c_str());
	for (int i = 0; i < files.size(); ++i)
	{
		cout << "Loading file: " << files[i].c_str() << endl;
		ifstream read_in(files[i].c_str());
		string temp;
		getline(read_in, temp);
		bool fq_flag = 0;
		if (temp[0] == '@')
		{
			fq_flag = 1;
		}
		read_in.close();

		
		string tag, n_tag, seq, qs;
		bool read_success = 1;
		read_in.open(files[i].c_str());
		
		while (read_success)
		{
			if (fq_flag)
			{
				read_success = get_a_fastq_read(read_in, tag, seq, qs);
			}
			else
			{
				read_success = get_a_fasta_read(read_in, tag, seq, n_tag);
			}
			if (read_success)
			{
				//n_reads++;
			}
			else
			{
				break;
			}
			
			if (seq.size() >= trim_len)
			{
				tag[0] = '>';
				o_selected_reads << tag << endl;
				o_selected_reads << seq << endl;
				sum += seq.size();
				n_reads++;
			}
		}
		
		read_in.close();
	}
	cout << n_reads << " reads selected." << endl;
	o_log << n_reads << " reads selected." << endl;
	cout << "Sum: "<<sum << endl;
	o_log << "Sum: " << sum << endl;
}