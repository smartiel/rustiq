mkdir -p results/chem
for file_name in $(ls chem/*.txt)
do
    ../target/release/rustiq $file_name --metric=count --onlyinfo >> ./results/$file_name;
    ../target/release/rustiq $file_name --metric=depth --onlyinfo >> ./results/$file_name;
    ../target/release/rustiq $file_name --metric=count --onlyinfo --keeporder >> ./results/$file_name;
    ../target/release/rustiq $file_name --metric=depth --onlyinfo --keeporder >> ./results/$file_name;
done