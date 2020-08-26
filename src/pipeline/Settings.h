#ifndef SETTINGS_H
#define SETTINGS_H

struct Settings {
    // Loads settings from disk, or fills in default values if no saved settings exist.
    void loadSettingsOrDefaults();



};

// The global Settings object, will be initialized by MainWindow
extern Settings settings;

#endif // SETTINGS_H
