#ifndef UTILS_HPP
#define UTILS_HPP

namespace Utils{

template <class NumType>
inline NumType StrToNum(const std::string& inStr, const NumType defaultValue = NumType(), bool* outSucceed = nullptr){
    std::istringstream iss (inStr);
    NumType result;
    iss >> result;
    if (iss.good() || iss.eof() ) {
        if(outSucceed){
            *outSucceed = true;
        }
        return result;
    }
    if(outSucceed){
        *outSucceed = false;
    }
    return defaultValue;
}

class ParamHelper{
    std::vector<std::string> params;
    bool hasFailed;

    int paramPosition(const std::vector<std::string>& inNames) const{
        for(int idxArg = 0 ; idxArg < int(params.size()) ; ++idxArg){
            for(auto& name : inNames){
                if(name == params[idxArg]){
                    return idxArg;
                }
            }
        }
        return -1;
    }

    int paramPosition(const std::initializer_list<std::string>& inNames) const{
        for(int idxArg = 0 ; idxArg < int(params.size()) ; ++idxArg){
            for(auto& name : inNames){
                if(name == params[idxArg]){
                    return idxArg;
                }
            }
        }
        return -1;
    }

public:
    ParamHelper(int argc, char** argv): hasFailed(false){
        params.insert(params.end(), argv, argv+argc);
    }

    bool parseHasFailed() const{
        return hasFailed;
    }

    int getNbParams() const{
        return int(params.size());
    }

    bool paramExist(std::initializer_list<std::string>& inNames) const{
        return paramPosition(inNames) != -1;
    }

    bool paramExist(const std::vector<std::string>& inNames) const{
        return paramPosition(inNames) != -1;
    }

    template <class ParamType>
    inline ParamType getValue(const std::vector<std::string>& inNames,
                       const ParamType inDefault) {
        const int idx = paramPosition(inNames);
        if(idx == -1 || idx == int(params.size()) - 1){
            if(idx == int(params.size()) - 1){
                hasFailed = true;
            }
            return inDefault;
        }
        bool succeed;
        ParamType value = StrToNum(params[idx+1], inDefault, &succeed);
        if(!succeed){
            hasFailed = true;
        }
        return value;
    }

    template <class ParamType>
    inline ParamType getValue(const std::initializer_list<std::string>& inNames,
                       const ParamType inDefault) {
        std::vector<std::string> names;
        names.insert(names.end(), inNames.begin(), inNames.end());
        return getValue<ParamType>(names, inDefault);
    }

    template <class ParamType>
    inline ParamType getValue(const std::vector<std::string>& inNames) {
        const int idx = paramPosition(inNames);
        if(idx == -1 || idx == int(params.size()) - 1){
            hasFailed = true;
            return ParamType();
        }
        bool succeed;
        ParamType value = StrToNum(params[idx+1], ParamType(), &succeed);
        if(!succeed){
            hasFailed = true;
        }
        return value;
    }

    template <class ParamType>
    inline ParamType getValue(const std::initializer_list<std::string>& inNames) {
        std::vector<std::string> names;
        names.insert(names.end(), inNames.begin(), inNames.end());
        return getValue<ParamType>(names);
    }

    inline std::string getStr(const std::vector<std::string>& inNames, const std::string inDefault) {
        const int idx = paramPosition(inNames);
        if(idx == -1 || idx == int(params.size()) - 1){
            if(idx == int(params.size()) - 1){
                hasFailed = true;
            }
            return inDefault;
        }
        return params[idx+1];
    }

    inline  std::string getStr(const std::initializer_list<std::string>& inNames, const std::string inDefault) {
        std::vector<std::string> names;
        names.insert(names.end(), inNames.begin(), inNames.end());
        return getStr(names, inDefault);
    }

    inline std::string getStr(const std::vector<std::string>& inNames) {
        const int idx = paramPosition(inNames);
        if(idx == -1 || idx == int(params.size()) - 1){
            hasFailed = true;
            return "";
        }
        return params[idx+1];
    }

    inline  std::string getStr(const std::initializer_list<std::string>& inNames) {
        std::vector<std::string> names;
        names.insert(names.end(), inNames.begin(), inNames.end());
        return getStr(names);
    }

    template <class ParamType>
    inline std::vector<ParamType> getValues(const std::vector<std::string>& inNames) {
        int idx = paramPosition(inNames);
        if(idx == -1 || idx == int(params.size()) - 1){
            hasFailed = true;
            return std::vector<ParamType>();
        }

        idx += 1;

        std::vector<ParamType> values;
        while(idx != int(params.size())){
            bool succeed;
            ParamType value = StrToNum(params[idx], ParamType(), &succeed);
            if(succeed == false){
                break;
            }
            values.push_back(value);
            idx += 1;
        }

        return values;
    }
};


}

#endif
