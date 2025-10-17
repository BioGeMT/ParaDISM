#!/bin/bash

# Interactive PKD1 Mapper Pipeline
# Provides a user-friendly interface for running the mapper with clean output

set -euo pipefail

# Colors and formatting
BLUE='\033[0;34m'
CYAN='\033[0;36m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BOLD='\033[1m'
DIM='\033[2m'
RESET='\033[0m'

# Print header
print_header() {
    clear
    echo -e "${BOLD}${BLUE}"
    cat << "EOF"
███╗   ███╗ █████╗ ██████╗ ██████╗ ███████╗██████╗
████╗ ████║██╔══██╗██╔══██╗██╔══██╗██╔════╝██╔══██╗
██╔████╔██║███████║██████╔╝██████╔╝█████╗  ██████╔╝
██║╚██╔╝██║██╔══██║██╔═══╝ ██╔═══╝ ██╔══╝  ██╔══██╗
██║ ╚═╝ ██║██║  ██║██║     ██║     ███████╗██║  ██║
╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝     ╚══════╝╚═╝  ╚═╝
EOF
    echo -e "${RESET}"
}

# Print section header
print_section() {
    echo "" >&2
    echo -e "${BOLD}${CYAN}▶ $1${RESET}" >&2
    echo "" >&2
}

# Smart file finder - looks for files with r1/r2/ref/sam in name
find_smart_file() {
    local search_term=$1  # "r1", "r2", "ref", or "sam"
    local extension=$2     # "*.fq", "*.fa", or "*.sam"

    # Case-insensitive search for files containing search term (current directory only)
    mapfile -t candidates < <(find . -maxdepth 1 -iname "*${search_term}*${extension#\*}" -type f 2>/dev/null | sort)

    if [ ${#candidates[@]} -eq 1 ]; then
        # Found exactly one match
        echo "${candidates[0]}"
        return 0
    else
        # No match or multiple matches
        return 1
    fi
}

# List files with numbers
list_files() {
    local pattern=$1
    local description=$2

    # Search current directory only
    mapfile -t files < <(find . -maxdepth 1 -name "$pattern" -type f 2>/dev/null | sort)

    if [ ${#files[@]} -eq 0 ]; then
        echo -e "${YELLOW}No $description found${RESET}" >&2
        return 1
    fi

    # Calculate width needed for numbering
    local max_num=${#files[@]}
    local num_width=${#max_num}

    echo "" >&2
    for i in "${!files[@]}"; do
        local num=$((i + 1))
        local file=${files[$i]}
        local size=$(du -h "$file" 2>/dev/null | cut -f1)
        printf "  ${CYAN}%${num_width}d${RESET} │ ${GREEN}%-50s${RESET} ${DIM}(%s)${RESET}\n" "$num" "$file" "$size" >&2
    done
    echo "" >&2
}

# Select file without smart detection (just show all matching pattern)
select_file_simple() {
    local pattern=$1       # "*.fa", "*.sam", etc.
    local description=$2   # User-friendly name

    # Show full list (current directory only)
    mapfile -t files < <(find . -maxdepth 1 -name "$pattern" -type f 2>/dev/null | sort)

    if [ ${#files[@]} -eq 0 ]; then
        echo -e "${RED}No $description found!${RESET}" >&2
        read -p "$(echo -e ${GREEN}Enter path manually: ${RESET})" custom_path </dev/tty
        echo "$custom_path"
        return 0
    fi

    # If exactly one file, auto-suggest it
    if [ ${#files[@]} -eq 1 ]; then
        local smart_file="${files[0]}"
        local size=$(du -h "$smart_file" 2>/dev/null | cut -f1)
        echo -e "${BOLD}Found:${RESET} ${GREEN}${smart_file}${RESET} ${DIM}(${size})${RESET}" >&2

        while true; do
            read -p "$(echo -e ${GREEN}Is this your $description? [Y/n]: ${RESET})" confirm </dev/tty
            case "$confirm" in
                ""|[Yy]|[Yy][Ee][Ss])
                    echo "$smart_file"
                    return 0
                    ;;
                [Nn]|[Nn][Oo])
                    echo -e "${RED}No other $description found!${RESET}" >&2
                    read -p "$(echo -e ${GREEN}Enter path manually: ${RESET})" custom_path </dev/tty
                    echo "$custom_path"
                    return 0
                    ;;
                *)
                    echo -e "${YELLOW}Please answer y or n${RESET}" >&2
                    ;;
            esac
        done
    fi

    # Multiple files - show numbered list
    list_files "$pattern" "$description"

    local max_num=${#files[@]}

    while true; do
        read -p "$(echo -e ${GREEN}Select 1-$max_num or ENTER for first: ${RESET})" choice </dev/tty

        if [ -z "$choice" ]; then
            choice=1
        fi

        if [[ "$choice" =~ ^[0-9]+$ ]] && [ "$choice" -ge 1 ] && [ "$choice" -le "$max_num" ]; then
            echo "${files[$((choice - 1))]}"
            return 0
        fi

        echo -e "${RED}Invalid selection. Enter 1-$max_num${RESET}" >&2
    done
}

# Select file with smart detection
select_file_smart() {
    local search_term=$1   # "r1", "r2", "ref", or "sam"
    local pattern=$2       # "*.fq", "*.fa", or "*.sam"
    local description=$3   # User-friendly name

    # Try smart detection first
    local smart_file
    if smart_file=$(find_smart_file "$search_term" "$pattern"); then
        # Found exactly one candidate
        local size=$(du -h "$smart_file" 2>/dev/null | cut -f1)
        echo -e "${BOLD}Found:${RESET} ${GREEN}${smart_file}${RESET} ${DIM}(${size})${RESET}" >&2

        while true; do
            read -p "$(echo -e ${GREEN}Is this your $description? [Y/n]: ${RESET})" confirm </dev/tty
            case "$confirm" in
                ""|[Yy]|[Yy][Ee][Ss])
                    echo "$smart_file"
                    return 0
                    ;;
                [Nn]|[Nn][Oo])
                    break
                    ;;
                *)
                    echo -e "${YELLOW}Please answer y or n${RESET}" >&2
                    ;;
            esac
        done
    fi

    # Smart detection failed or user said no - show full list (current directory only)
    mapfile -t files < <(find . -maxdepth 1 -name "$pattern" -type f 2>/dev/null | sort)

    if [ ${#files[@]} -eq 0 ]; then
        echo -e "${RED}No $description found!${RESET}" >&2
        read -p "$(echo -e ${GREEN}Enter path manually: ${RESET})" custom_path </dev/tty
        echo "$custom_path"
        return 0
    fi

    list_files "$pattern" "$description"

    local max_num=${#files[@]}

    while true; do
        if [ $max_num -eq 1 ]; then
            read -p "$(echo -e ${GREEN}Press ENTER to select: ${RESET})" choice </dev/tty
        else
            read -p "$(echo -e ${GREEN}Select 1-$max_num or ENTER for first: ${RESET})" choice </dev/tty
        fi

        if [ -z "$choice" ]; then
            choice=1
        fi

        if [[ "$choice" =~ ^[0-9]+$ ]] && [ "$choice" -ge 1 ] && [ "$choice" -le "$max_num" ]; then
            echo "${files[$((choice - 1))]}"
            return 0
        fi

        echo -e "${RED}Invalid selection. Enter 1-$max_num${RESET}" >&2
    done
}

# Select aligner
select_aligner() {
    print_section "Select Read Aligner"

    echo -e "  ${CYAN} 1${RESET} │ ${GREEN}BWA-MEM2${RESET}     ${DIM}(Default, fast for short reads)${RESET}" >&2
    echo -e "  ${CYAN} 2${RESET} │ ${GREEN}Bowtie2${RESET}      ${DIM}(Alternative for short reads)${RESET}" >&2
    echo -e "  ${CYAN} 3${RESET} │ ${GREEN}Minimap2${RESET}     ${DIM}(Versatile, supports long reads)${RESET}" >&2
    echo "" >&2

    while true; do
        read -p "$(echo -e ${GREEN}Select aligner [1-3] or press ENTER for BWA-MEM2: ${RESET})" choice </dev/tty

        case "$choice" in
            ""|1)
                echo "bwa-mem2"
                return 0
                ;;
            2)
                echo "bowtie2"
                return 0
                ;;
            3)
                # Ask for minimap2 profile
                print_section "Select Minimap2 Profile"
                echo -e "  ${CYAN} 1${RESET} │ ${GREEN}short${RESET}     ${DIM}(Illumina paired-end)${RESET}" >&2
                echo -e "  ${CYAN} 2${RESET} │ ${GREEN}pacbio${RESET}    ${DIM}(PacBio HiFi/CLR)${RESET}" >&2
                echo -e "  ${CYAN} 3${RESET} │ ${GREEN}nanopore${RESET}  ${DIM}(Oxford Nanopore)${RESET}" >&2
                echo "" >&2

                while true; do
                    read -p "$(echo -e ${GREEN}Select profile [1-3] or press ENTER for short: ${RESET})" profile_choice </dev/tty

                    case "$profile_choice" in
                        ""|1)
                            MINIMAP2_PROFILE="short"
                            echo "minimap2"
                            return 0
                            ;;
                        2)
                            MINIMAP2_PROFILE="pacbio"
                            echo "minimap2"
                            return 0
                            ;;
                        3)
                            MINIMAP2_PROFILE="nanopore"
                            echo "minimap2"
                            return 0
                            ;;
                        *)
                            echo -e "${RED}Invalid selection. Please enter 1-3${RESET}" >&2
                            ;;
                    esac
                done
                ;;
            *)
                echo -e "${RED}Invalid selection. Please enter 1-3${RESET}" >&2
                ;;
        esac
    done
}

# Select thread count
select_threads() {
    local default_threads=4
    local max_threads=$(nproc)

    read -p "$(echo -e ${GREEN}Number of threads [1-$max_threads] or press ENTER for $default_threads: ${RESET})" threads </dev/tty

    if [ -z "$threads" ]; then
        echo $default_threads
        return 0
    fi

    if [[ "$threads" =~ ^[0-9]+$ ]] && [ "$threads" -ge 1 ] && [ "$threads" -le "$max_threads" ]; then
        echo $threads
        return 0
    else
        echo -e "${YELLOW}Invalid thread count. Using default: $default_threads${RESET}" >&2
        echo $default_threads
        return 0
    fi
}

# Ask if user wants to use existing SAM file
ask_use_sam() {
    echo "" >&2

    while true; do
        read -p "$(echo -e ${GREEN}Do you have an existing SAM alignment file? [y/N]: ${RESET})" choice </dev/tty

        case "$choice" in
            ""|[Nn]|[Nn][Oo])
                echo ""
                return 0
                ;;
            [Yy]|[Yy][Ee][Ss])
                echo "sam"
                return 0
                ;;
            *)
                echo -e "${YELLOW}Please answer y or n${RESET}" >&2
                ;;
        esac
    done
}

# Main script
main() {
    print_header

    # Select R1
    print_section "Select Forward Reads (R1)"
    R1=$(select_file_smart "r1" "*.fq" "R1 file")
    [ -z "$R1" ] && exit 1
    echo -e "${GREEN}✓${RESET} Selected R1: ${BOLD}$R1${RESET}"

    # Select R2
    print_section "Select Reverse Reads (R2)"
    R2=$(select_file_smart "r2" "*.fq" "R2 file")
    [ -z "$R2" ] && exit 1
    echo -e "${GREEN}✓${RESET} Selected R2: ${BOLD}$R2${RESET}"

    # Select Reference
    print_section "Select Reference"
    REF=$(select_file_simple "*.fa" "reference file")
    [ -z "$REF" ] && exit 1
    echo -e "${GREEN}✓${RESET} Selected Reference: ${BOLD}$REF${RESET}"

    # Ask about SAM file
    USE_SAM=$(ask_use_sam)

    if [ "$USE_SAM" = "sam" ]; then
        # Select existing SAM file
        print_section "Select SAM File"
        SAM=$(select_file_smart "sam" "*.sam" "SAM file")
        [ -z "$SAM" ] && exit 1
        echo -e "${GREEN}✓${RESET} Selected SAM: ${BOLD}$SAM${RESET}"
    else
        # Select Aligner (print_section is inside select_aligner)
        ALIGNER=$(select_aligner)
        if [ "$ALIGNER" = "minimap2" ]; then
            echo -e "${GREEN}✓${RESET} Selected Aligner: ${BOLD}$ALIGNER ($MINIMAP2_PROFILE)${RESET}"
        else
            echo -e "${GREEN}✓${RESET} Selected Aligner: ${BOLD}$ALIGNER${RESET}"
        fi

        # Select Threads
        print_section "Configure Resources"
        THREADS=$(select_threads)
        echo -e "${GREEN}✓${RESET} Using threads: ${BOLD}$THREADS${RESET}"
    fi

    # Confirm and run
    echo ""
    echo -e "${BOLD}${YELLOW}═══════════════════════════════════════════════════════${RESET}"
    echo -e "${BOLD}Ready to start pipeline${RESET}"
    echo -e "${YELLOW}═══════════════════════════════════════════════════════${RESET}"
    echo ""

    read -p "$(echo -e ${GREEN}Press ENTER to start or Ctrl+C to cancel...${RESET})"

    # Run pipeline
    print_section "Running Pipeline"

    # Build mapper.sh arguments
    MAPPER_ARGS="--read1 $R1 --read2 $R2 --reference $REF"

    if [ "$USE_SAM" = "sam" ]; then
        MAPPER_ARGS="$MAPPER_ARGS --sam $SAM"
    else
        MAPPER_ARGS="$MAPPER_ARGS --aligner $ALIGNER --threads $THREADS"
        if [ "$ALIGNER" = "minimap2" ]; then
            MAPPER_ARGS="$MAPPER_ARGS --minimap2-profile $MINIMAP2_PROFILE"
        fi
    fi

    echo -e "${DIM}Output directory: ./output${RESET}"
    echo ""

    # Run the pipeline (it handles MSA and everything else)
    echo -e "${CYAN}Running mapper pipeline...${RESET}"
    echo ""
    MAPPER_CALLED_FROM_INTERACTIVE=1 ./mapper.sh $MAPPER_ARGS

    # Success
    echo ""
    echo -e "${BOLD}${GREEN}═══════════════════════════════════════════════════════${RESET}"
    echo -e "${BOLD}${GREEN}Pipeline Complete!${RESET}"
    echo -e "${GREEN}═══════════════════════════════════════════════════════${RESET}"
    echo ""
    echo -e "${BOLD}Outputs:${RESET}"
    echo -e "  ${CYAN}•${RESET} Unique mappings: ${GREEN}./output/unique_mappings.tsv${RESET}"
    echo -e "  ${CYAN}•${RESET} Gene-specific FASTQs: ${GREEN}./output/fastq/${RESET}"
    echo -e "  ${CYAN}•${RESET} Gene-specific BAMs: ${GREEN}./output/bam/${RESET}"
    echo ""
}

# Run main
main
