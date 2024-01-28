from transformers import BertPreTrainedModel, BertEncoder
from torch import nn
import torch

class BertClassifier(BertPreTrainedModel):

    def __init__(self, config, add_pooling_layer=True):
        super().__init__(config)
        self.config = config

        self.encoder = BertEncoder(config)

        self.pooler = BertPooler(config) if add_pooling_layer else None

        classifier_dropout = (
            config.classifier_dropout if config.classifier_dropout is not None else config.hidden_dropout_prob
        )
        self.dropout = nn.Dropout(classifier_dropout)
        self.classifier = nn.Linear(config.hidden_size, config.num_labels)

        self.post_init()
    
    def forward(
        self,
        input_ids,
        attention_mask,
    ):
        
        input_shape = input_ids.size()

        batch_size, seq_length = input_shape

        extended_attention_mask: torch.Tensor = self.get_extended_attention_mask(attention_mask, input_shape)

        encoder_outputs = self.encoder(
            input_ids,
            attention_mask=extended_attention_mask
        )

        sequence_output = encoder_outputs[0]
        pooled_output = self.pooler(sequence_output) if self.pooler is not None else None

        pooled_output = self.dropout(pooled_output)
        logits = self.classifier(pooled_output)

        return logits