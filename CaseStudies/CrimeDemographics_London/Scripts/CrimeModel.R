sub = crime[,c(-2:-4,-7)]
train.crime = subset(sub, year < 2016)
test.crime = subset(sub, year >= 2016)

# Data Reduction
#--------------------------
lsoa_enc = as.numeric(factor(train.crime$lsoa_code))
train.enc = cbind(lsoa_enc,train.crime[,2:3])
train.melt = melt(train.enc, id=c("lsoa_enc","year"))
train.years = dcast(train.melt, lsoa_enc+year~variable, sum)

lsoa_enc2 = as.numeric(factor(test.crime$lsoa_code))
test.enc = cbind(lsoa_enc2,test.crime[,2:3])
test.melt = melt(test.enc, id=c("lsoa_enc2","year"))
test.years = dcast(test.melt, lsoa_enc2+year~variable, sum)

# LSOA Encoding Key
#-------------------------
enc.bind = cbind(data.frame(train.crime$lsoa_code), data.frame(lsoa_enc))
key = melt(enc.bind, id = c("train.crime.lsoa_code", "lsoa_enc"))
primary.key = dcast(key, train.crime.lsoa_code + lsoa_enc ~ .)
rm(key,enc.bind,lsoa_enc)
primary.key = subset(primary.key[,-3])