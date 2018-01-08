#include "include/showusage.h"
#include <string>
#include "include/config.h"

void show_usage(int argc, char* argv[], char* config_file) {
    bool valid_args = false;
    do {
        Config* pcfg = Config::get_instance();
        if(!pcfg->init(config_file))
              break;         

        if (argc < 3) break;
        int param_idx = 1;

        std::string up_down_flag = pcfg->get_up_down_flag();
        std::string file_name = pcfg->get_file_name();
        std::string dst_path = pcfg->get_dst_path();
        std::string ftp_host = pcfg->get_ftp_host();
        bool is_rename = pcfg->get_is_rename();
        std::string new_filename = pcfg->get_new_filename();

        if (argc > param_idx) {
            up_down_flag = std::string(argv[param_idx]);
            pcfg->set_up_down_flag(up_down_flag);
        }
        param_idx++;

        if (argc > param_idx) {
            file_name = std::string(argv[param_idx]);
            pcfg->set_file_name(file_name);
        }
        param_idx++;

        if (argc > param_idx) {
            ftp_host = std::string(argv[param_idx]);
            pcfg->set_ftp_host(ftp_host);
        }
            
        if(up_down_flag == "download")
        {
           param_idx++;

           if (argc > param_idx) {
               dst_path = std::string(argv[param_idx]);
               pcfg->set_dst_path(dst_path);
           }
        }
        param_idx++;
        if (argc > param_idx) {
            is_rename = bool(argv[param_idx]);
            pcfg->set_is_rename(is_rename);
           
            param_idx++;
            if (is_rename == 1 && argc > param_idx) {
               new_filename = std::string(argv[param_idx]);
               pcfg->set_new_filename(new_filename);
            }
        }
        else
        {
            is_rename = 0;
            pcfg->set_is_rename(is_rename);

        }


        valid_args = true;

    } while (0);

    if (!valid_args) {
        printf("usage:\n");
        printf("./seeta_ftp up_down_flag file_name ftp_host dst_path(optional) is_rename(optional) new_filename(optional)\n");
        exit(EXIT_SUCCESS);
    }
}
