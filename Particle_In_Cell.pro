SOURCES += \
    window.cpp \
    main.cpp \
    plasma.cpp

HEADERS += \
    window.h \
    plasma.h

QT += \
    core gui \
    widgets \
    charts

INCLUDEPATH += \
    /opt/local/include

LIBS += \
    -L /opt/local/lib \
    -lfftw3 \
    -lm
