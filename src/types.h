#pragma once
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

enum class BlockType { EDGE, HARD_MACRO, SOFT };

struct FTConversion {
    double rate[4];
};

struct Block {
    std::string name;
    BlockType type;
    double area;          
    double width, height; 
    double min_ar, max_ar;
    std::vector<std::string> locations; 
    FTConversion ft;
    bool has_fixed_wh;    

    double lx = 0, ly = 0;

    double get_ar() const { return (height > 0) ? width / height : 1.0; }

    // Get the target area after considering FT expansion
    double get_target_area(int ft_nets) const {
        if (ft_nets <= 0) return area;
        double rate = ft_conversion_rate(ft_nets);
        double base_side = std::sqrt(area);
        double extend = ((double)ft_nets / 25.0) * rate / 2.0;
        double new_side = base_side + extend;
        return new_side * new_side;
    }

    double ft_conversion_rate(int ft_nets) const {
        if (ft_nets <= 3000)       return ft.rate[0];
        else if (ft_nets <= 6000)  return ft.rate[1];
        else if (ft_nets <= 9000)  return ft.rate[2];
        else                       return ft.rate[3];
    }
};

struct Channel {
    std::string name;
    double lx, ly, width, height;
    
    int net_count = 0; 
    double nets_x = 0; 
    double nets_y = 0;

    double cap_x() const { return height * 25.0; } 
    double cap_y() const { return width * 25.0; }  

    bool overflowed() const { return nets_x > cap_x() || nets_y > cap_y(); }
};

struct Connection {
    int from, to;   
    int nets;
};

struct PathSegment {
    std::string rect_name; 
    int edge_in;  
    int edge_out;
};

struct RoutePath {
    int nets;
    std::vector<PathSegment> segments;
};

struct Outline {
    double max_width, max_height;
    double cur_width, cur_height; 
};

struct Design {
    Outline outline;
    std::vector<Block> blocks;
    std::vector<Connection> connections; 
    std::vector<Channel> channels;
    std::vector<RoutePath> paths;
    double alpha = 1.0;

    int block_idx(const std::string& name) const {
        for (int i = 0; i < (int)blocks.size(); i++)
            if (blocks[i].name == name) return i;
        return -1;
    }
};