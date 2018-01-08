// Copyright 2017-2018 SeetaTech
// Editor BoxuanXu
// Time 2017/12/6
#include "include/config.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#ifdef WITH_FILECONF
#include <json/json.h>
#endif
//#include "include/eupulogger4system.h"

namespace {
static const char* FTP_PATH = "ftp_path";
static const char* USERNAME = "username";
static const char* PASSWORD = "password";
}  // namespace

Config::Config() {}

Config::~Config() {}

bool Config::init(const std::string& cfg_path) {
    bool ret = false;
    do {
#ifdef WITH_FILECONF
        std::ifstream f;
        f.open(cfg_path.c_str());
        if (!f.is_open()) {
            printf("Config::init() open config file %s failed",
                cfg_path.c_str());
            break;
        }

        Json::Reader r;
        Json::Value v;
        if (!r.parse(f, v, false)) {
            printf("Config::init() parse config %s failed",
                cfg_path.c_str());
            break;
        }

        cfg_.ftp_path = v[FTP_PATH].asString();
        cfg_.username = v[USERNAME].asString();
        cfg_.password = v[PASSWORD].asString();
#else
        cfg_.ftp_path = "guest_dir/";
        cfg_.username = "guest";
        cfg_.password = "guest";
        
#endif
        cfg_.up_down_flag = "upload";
        cfg_.file_name = "";
        cfg_.dst_path = ".";
        cfg_.ftp_host = "192.168.1.247";
        cfg_.is_rename = 0;
        cfg_.new_filename = "";  
        ret = true;
    } while (0);
    return ret;
}

void Config::output() {
    if(cfg_.up_down_flag == "upload")
    {
        printf("up_down_flag: %s\n", cfg_.up_down_flag.c_str());
        printf("src_file: %s\n", cfg_.file_name.c_str());
        printf("ftp_host: %s\n", cfg_.ftp_host.c_str());
        printf("%s: %s\n", FTP_PATH,cfg_.ftp_path.c_str());
        printf("%s: %s\n", USERNAME,cfg_.username.c_str());
        printf("%s: %s\n", PASSWORD,cfg_.password.c_str());
        printf("is_rename: %d\n", cfg_.is_rename);
        printf("new_filename: %s\n", cfg_.new_filename.c_str());
    }
    else
    {
        printf("up_down_flag: %s\n", cfg_.up_down_flag.c_str());
        printf("src_file: %s\n", cfg_.file_name.c_str());
        printf("dst_path: %s\n", cfg_.dst_path.c_str());
        printf("ftp_host: %s\n", cfg_.ftp_host.c_str());
        printf("%s: %s\n", FTP_PATH,cfg_.ftp_path.c_str());
        printf("%s: %s\n", USERNAME,cfg_.username.c_str());
        printf("%s: %s\n", PASSWORD,cfg_.password.c_str());
        printf("is_rename: %d\n", cfg_.is_rename);
        printf("new_filename: %s\n", cfg_.new_filename.c_str());

    } 
}
