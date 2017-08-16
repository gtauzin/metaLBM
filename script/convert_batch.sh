ls energy-removal_energy-spectrum_* | awk 'BEGIN{ a=0 }{ printf "convert %s %d.png\n", $0, a++ }' | bash
