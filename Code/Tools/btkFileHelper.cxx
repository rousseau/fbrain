#include "btkFileHelper.h"


namespace btk
{
//------------------------------------------------------------------------------------------------
bool FileHelper::FileExist(std::string &file)
{
    std::fstream fs;

    fs.open(file.c_str());

    return fs.is_open();
}
//------------------------------------------------------------------------------------------------
bool FileHelper::FileExistPrintIt(std::string &file)
{
    std::fstream fs;

    fs.open(file.c_str());

    if(fs.is_open() == true)
    {
        std::cout<<"File "<<file<<" exist."<<std::endl;
        return true;
    }
    else
    {
        std::cout<<"File "<<file<<" does not exist."<<std::endl;
        return false;
    }
}
//------------------------------------------------------------------------------------------------
void FileHelper::FilesExist(std::vector<std::string> &files, std::vector<bool> &result)
{
    result.resize(files.size());
    std::fstream fs;

    for(int i = 0; i< files.size(); i++)
    {
        fs.open(files[i].c_str());

        if(fs.is_open() == true)
        {
            result[i] = true;
        }
        else
        {
            result[i] = false;
        }

    }
}
//------------------------------------------------------------------------------------------------
void FileHelper::FilesExistPrintIt(std::vector<std::string> &files, std::vector<bool> &result)
{
    result.resize(files.size());
    std::fstream fs;

    for(int i = 0; i< files.size(); i++)
    {
        fs.open(files[i].c_str());

        if(fs.is_open() == true)
        {
            result[i] = true;
            std::cout<<"File "<<files[i]<<" exist."<<std::endl;
        }
        else
        {
            result[i] = false;
            std::cout<<"File "<<files[i]<<" does not exist."<<std::endl;
        }

    }
}
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
}
