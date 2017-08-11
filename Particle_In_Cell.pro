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


LIBS += \
    -lfftw3 \
    -lm
