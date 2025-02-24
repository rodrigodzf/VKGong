/*
 * VK-Gong code 
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 */

#include "SettingsManager.h"
#include "GlobalSettings.h"
#include "Logger.h"

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cctype>
using namespace std;

SettingsManager *SettingsManager::instance = NULL;

SettingsManager::SettingsManager()
{
}

SettingsManager::~SettingsManager()
{
}


void SettingsManager::putSettingInt(string compname, string settingname, int value)
{
    char buffer[100];
    sprintf(buffer, "%d", value);
    putSetting(compname, settingname, buffer);
}

void SettingsManager::putSetting(string compname, string settingname, double value)
{
    char buffer[100];
    sprintf(buffer, "%.20f", value);
    putSetting(compname, settingname, buffer);
}

void SettingsManager::putSettingBool(string compname, string settingname, bool value)
{
    if (value) putSetting(compname, settingname, "true");
    else putSetting(compname, settingname, "false");
}

void SettingsManager::putSetting(string compname, string settingname, string value)
{
    Setting setting;
    setting.component = compname;
    setting.setting = settingname;
    setting.value = value;
    settings.push_back(setting);
    logMessage(1, "Adding local setting %s for component %s, value %s",
	       setting.setting.c_str(), setting.component.c_str(),
	       setting.value.c_str());
}

bool SettingsManager::putSetting(char *arg1, char *arg2)
{
    bool result = false;
    string str = arg1;
    int colon = str.find(":");
    string settingname = str.substr(1, colon-1);
    string compname = str.substr(colon+1);

    string value;
    if ((arg2 == NULL) ||
	((arg2[0] == '-') && (!isdigit(arg2[1])))) {
	// next arg either doesn't exist or is start of a different option
	value = "true";
    }
    else {
	value = arg2;
	result = true;
    }
    
    size_t comma = compname.find(",");
    if (comma != string::npos) {
	// handle multiple component case
	size_t start = 0;
	while (comma != string::npos) {
	    string thiscomp = compname.substr(start, comma - start);
	    putSetting(thiscomp, settingname, value);
	    start = comma + 1;
	    comma = compname.find(",", start);
	}
	string finalcomp = compname.substr(start);
	putSetting(finalcomp, settingname, value);
	return result;
    }


    Setting setting;
    setting.component = compname;
    setting.setting = settingname;
    setting.value = value;
    settings.push_back(setting);

    logMessage(1, "Adding local setting %s for component %s, value %s",
	       setting.setting.c_str(), setting.component.c_str(),
	       setting.value.c_str());
    return result;
}

Setting *SettingsManager::findSetting(string component, string setting)
{
    int i;
    for (i = 0; i < settings.size(); i++) {
	Setting *s = &(settings.at(i));

	if ((s->component == component) && (s->setting == setting)) {
	    return s;
	}
    }
    return NULL;
}

string SettingsManager::getStringSetting(string component, string setting)
{
    Setting *s = findSetting(component, setting);
    if (!s) return "";
    return s->value;
}

double SettingsManager::getDoubleSetting(string component, string setting)
{
    Setting *s = findSetting(component, setting);
    if (!s) {
	GlobalSettings *gs = GlobalSettings::getInstance();
    // parse double
    }
    return atof(s->value.c_str());
}

int SettingsManager::getIntSetting(string component, string setting)
{
    Setting *s = findSetting(component, setting);
    if (!s) {
	GlobalSettings *gs = GlobalSettings::getInstance();

	return 0;
    }
    // parse int
    return atoi(s->value.c_str());
}

bool SettingsManager::getBoolSetting(string component, string setting)
{
    Setting *s = findSetting(component, setting);
    if (!s) {
	GlobalSettings *gs = GlobalSettings::getInstance();
	
	return false;
    }
    // parse boolean
    return (s->value == "true");
}
