#!/bin/bash

# Stoppt das Skript sofort, wenn ein Befehl fehlschlägt
set -e

# --- Standardwerte für die Flags ---
ACADOS_WITH_OPENMP="OFF"
ACADOS_WITH_OSQP="OFF"
ACADOS_WITH_QPOASES="OFF"
ACADOS_PYTHON="OFF"
CREATE_VENV="OFF"
AUTO_EXPORT="OFF"
PIP_EDITABLE="OFF"
VENV_NAME=".acados_env" # Standard-Venv-Name


# --- Hilfefunktion ---
show_help() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  --omp                     Aktiviert die OpenMP-Unterstützung."
    echo "  --osqp                    Aktiviert den OSQP-Solver."
    echo "  --qpoases                 Aktiviert den qpOASES-Solver."
    echo "  -p, --python              Aktiviert die Python-Schnittstelle."
    echo "  -v, --venv                Erstellt eine Python venv für die Installation (benötigt -p)."
    echo "  --venv-name <NAME>        Gibt den Namen für die venv an (Standard: .acados_env)."
    echo "  -e, --export              Fügt die nötigen export-Befehle automatisch zu Ihrer Shell-Konfigurationsdatei hinzu."
    echo "  -pe, --pip-editable       Macht das acados editable in pip."
    echo "  -h, --help                Zeigt diese Hilfe an."
    echo ""
    echo "Beispiel: $0 --qpoases -p -v --venv-name my_acados_env -e"
}

# --- Abhängigkeiten prüfen und installieren ---
check_and_install_dependencies() {
    echo "Überprüfe und installiere Abhängigkeiten..."
    command_exists() {
        command -v "$1" &> /dev/null
    }

    if command_exists apt-get; then
        REQUIRED_PACKAGES="git cmake build-essential python3-dev"
        if [[ "$CREATE_VENV" == "ON" ]]; then
            REQUIRED_PACKAGES="$REQUIRED_PACKAGES python3-venv"
        fi

        echo "Benötigte Pakete für Debian/Ubuntu: $REQUIRED_PACKAGES"
        echo "Aktualisiere Paketlisten (sudo erforderlich)..."
        sudo apt-get update -qq
        echo "Installiere Pakete..."
        sudo apt-get install -y $REQUIRED_PACKAGES

    elif command_exists brew; then
        echo "Benötigte Pakete für macOS: git, cmake, python"
        brew install git cmake python

    else
        echo "WARNUNG: Konnte keinen unterstützten Paketmanager (apt, brew) finden."
        echo "   Bitte stelle sicher, dass die folgenden Abhängigkeiten manuell installiert sind:"
        echo "   - git"
        echo "   - cmake"
        echo "   - Ein C/C++ Compiler (z.B. build-essential oder Xcode Command Line Tools)"
        echo "   - python3-dev"
        if [[ "$CREATE_VENV" == "ON" ]]; then
            echo "   - python3-venv"
        fi
    fi
    echo "Abhängigkeiten sind bereit."
}

# --- Argumente parsen ---
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --omp) ACADOS_WITH_OPENMP="ON";;
        --osqp) ACADOS_WITH_OSQP="ON";;
        --qpoases) ACADOS_WITH_QPOASES="ON";;
        -p|--python) ACADOS_PYTHON="ON";;
        -v|--venv) CREATE_VENV="ON";;
        --venv-name)
            if [[ -n "$2" && "$2" != --* ]]; then
                VENV_NAME="$2"
                shift
            else
                echo "Fehler: --venv-name benötigt einen Namen als Argument." >&2
                exit 1
            fi
            ;;
        -e|--export) AUTO_EXPORT="ON";;
        -pe|--pip-editable) PIP_EDITABLE="ON";;
        -h|--help) show_help; exit 0;;
        *) echo "Unbekannte Option: $1"; show_help; exit 1;;
    esac
    shift # Geht zur nächsten Option
done

# --- Logik-Prüfung ---
if [[ "$CREATE_VENV" == "ON" && "$ACADOS_PYTHON" == "OFF" ]]; then
    echo "Fehler: Die --venv Option kann nur zusammen mit -p oder --python verwendet werden."
    exit 1
fi

# --- Abhängigkeiten installieren (NEU) ---
check_and_install_dependencies

# --- Installationsprozess ---
ACADOS_INSTALL_DIR="$( cd "$( dirname "$0" )" && pwd )"

echo "Starte die Installation von acados..."

if [ -d "$ACADOS_INSTALL_DIR/acados" ]; then
    read -p "Das Verzeichnis 'acados' existiert bereits. Löschen und neu beginnen? (j/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Jj]$ ]]; then
        rm -rf "$ACADOS_INSTALL_DIR/acados"
    else
        echo "Installation abgebrochen."
        exit 1
    fi
fi

# 1. Acados klonen
git submodule update --init --recursive

# 2. Build und Kompilierung
mkdir -p build
cd build

echo "Konfiguriere das Build-System mit CMake..."
cmake .. \
    -DACADOS_WITH_OPENMP=${ACADOS_WITH_OPENMP} \
    -DACADOS_WITH_OSQP=${ACADOS_WITH_OSQP} \
    -DACADOS_WITH_QPOASES=${ACADOS_WITH_QPOASES} \
    -DACADOS_PYTHON=${ACADOS_PYTHON}

CPU_CORES=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)
echo "Kompiliere und installiere acados mit $CPU_CORES Kernen..."
make install -j$CPU_CORES

# 3. Python-Interface installieren
if [[ "$ACADOS_PYTHON" == "ON" ]]; then
    if [[ "$CREATE_VENV" == "ON" ]];
    then
        VENV_PATH="$HOME/$VENV_NAME"
        echo "Erstelle Python virtual environment in: $VENV_PATH"
        python3 -m venv "$VENV_PATH"
        source "$VENV_PATH/bin/activate"
    fi

    echo "Installiere Python-Schnittstelle (acados_template)..."
    cd $ACADOS_INSTALL_DIR/interfaces/acados_template

    if [[ "$CREATE_VENV" == "ON" ]]; then
        if [[ "$PIP_EDITABLE" == "ON" ]]; then
            pip install -e .
        else
            pip install .
        fi
    else
        echo "WARNUNG: Installation erfolgt systemweit. Nutze --venv für eine isolierte Installation."
        sudo pip install .
    fi
fi

# --- Abschließende Anweisungen ---
echo ""
echo "Acados-Installation abgeschlossen!"

if [[ "$AUTO_EXPORT" == "ON" ]]; then
    RC_FILE=""
    if [[ "$SHELL" == *"bash"* ]]; then
        RC_FILE="$HOME/.bashrc"
    elif [[ "$SHELL" == *"zsh"* ]]; then
        RC_FILE="$HOME/.zshrc"
    fi

    if [ -n "$RC_FILE" ]; then
        if ! grep -q "ACADOS_SOURCE_DIR" "$RC_FILE"; then
            echo "Füge acados-Pfade zu $RC_FILE hinzu..."
            {
                echo ""
                echo "# --- acados environment variables (hinzugefügt durch Installationsskript) ---"
                echo "export ACADOS_SOURCE_DIR=\"$ACADOS_INSTALL_DIR\""
                echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$ACADOS_SOURCE_DIR/lib"
            } >> "$RC_FILE"
            echo "   -> Erfolgreich!"
        else
            echo "☑️  acados-Pfade sind bereits in $RC_FILE vorhanden."
        fi
        echo "Bitte starte eine neue Terminalsitzung oder führe 'source $RC_FILE' aus, um die Änderungen zu aktivieren."
    else
        echo "Unbekannte Shell: $SHELL. Konnte Pfade nicht automatisch exportieren."
        echo "Bitte füge die folgenden Zeilen manuell zu deiner Shell-Konfigurationsdatei hinzu:"
        echo "export ACADOS_SOURCE_DIR=\"$ACADOS_INSTALL_DIR\""
        echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$ACADOS_SOURCE_DIR/lib"
    fi
else
    echo "Bitte füge die folgenden Zeilen manuell zu deiner ~/.bashrc oder ~/.zshrc hinzu:"
    echo "export ACADOS_SOURCE_DIR=\"$ACADOS_INSTALL_DIR\""
    echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$ACADOS_SOURCE_DIR/lib"
fi

if [[ "$CREATE_VENV" == "ON" ]]; then
    echo ""
    echo "Die Python-Pakete wurden in der virtuellen Umgebung '$VENV_NAME' installiert."
    echo "   Aktiviere sie vor der Verwendung mit:"
    echo "   source \"$HOME/$VENV_NAME/bin/activate\""
fi